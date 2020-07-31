#!/usr/bin/env python

"""Primary model code for hydronetworks."""

from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Point, shape
import rasterio
import geopandas as gpd
from pyproj import CRS
from pyproj.proj import Proj

from .streamorder import shreve, strahler
from .runoff import rainfall_runoff, discharge

wgs84 = CRS("EPSG:4326")
merc = CRS("EPSG:3395")
merc_pro = Proj("EPSG:3395")


class HydroNetworks:
    def __init__(self, rivers_path, dem_path, flow_path, flowacc_path):
        self.rivers = gpd.read_file(rivers_path)
        # project rivers to Mercator for distances and integer positioning
        self.rivers = self.rivers.to_crs(merc)
        self.dem = rasterio.open(dem_path)
        self.flow = rasterio.open(flow_path)
        self.flowacc = rasterio.open(flowacc_path)

    def save_results(self, path):
        if isinstance(path, str):
            path = Path(path)
        path.mkdir(exist_ok=True, parents=True)

        rivers_path = path / "rivers_out.geojson"
        nodes_path = path / "nodes_out.geojson"
        hydro_path = path / "hydro_points.geojson"

        self.rivers_out.to_file(rivers_path, driver="GeoJSON")
        self.nodes_out.drop("arcs", axis=1).to_file(nodes_path, driver="GeoJSON")
        self.hydro_points.to_file(hydro_path, driver="GeoJSON")

    def create_network(self):
        """
        Create the network and nodes from the river layer, ready to be stream ordered

        Structure for network:
        0   index
        1   fx
        2   fy
        3   lx
        4   ly
        5   node index first point
        6   node index last point
        7   stream order
        8   arc length

        Structure for nodes:
        0   index
        1   x
        2   y
        3   {} metadata
        4..   arc indices...
        """

        # import the start and end end points into the network array of arcs
        # and simultaneously import them into a separate node array (two nodes per arc)
        # there will be duplicate nodes, but they duplicates are ignored later
        counter = 0
        length = len(self.rivers)
        self.network = np.empty([length, 9], dtype=np.int32)
        self.nodes = []

        for index, row in self.rivers.iterrows():
            # the geometry variables (int so that they actual match when comparing)
            fx = int(row.geometry[0].xy[0][0])
            fy = int(row.geometry[0].xy[1][0])
            lx = int(row.geometry[0].xy[0][-1])
            ly = int(row.geometry[0].xy[1][-1])

            # add arc_length as a determinant of how much water contributes downstream
            # not sure if it's a good idea going forward...
            arc_length = int(row.geometry.length)

            # store the index as a column for easier access
            # the last column is for stream order
            self.network[counter] = [counter, fx, fy, lx, ly, -99, -99, -99, arc_length]

            # create nodes for the start point and end point separately
            self.nodes.append([2 * counter, fx, fy, {}])
            self.nodes.append([2 * counter + 1, lx, ly, {}])

            counter += 1

        # the most time consuming part of the script
        # runs through every arc and every node and where there is a match, add
        #  a reference to the index into the network array
        # once both nodes have been found for an arc, break that for loop and
        #  start with the next arc
        for arc in self.network:
            match_f, match_l = False, False
            for node in self.nodes:
                # if the fx and fy values match
                if arc[1] == node[1] and arc[2] == node[2]:
                    # add an index
                    arc[5] = node[0]
                    match_f = True
                # if the lx and ly values match
                elif arc[3] == node[1] and arc[4] == node[2]:
                    # add an index
                    arc[6] = node[0]
                    match_l = True
                if match_f and match_l > 1:
                    # reset and skip to the next arc
                    break

        # for every node, add references to every arc that connects to it
        for arc in self.network:
            # tell the arc's starting node that it exists
            self.nodes[arc[5]].append(arc[0])
            # and tell the arc's ending node that it exists
            self.nodes[arc[6]].append(arc[0])

    def assign_streamorder(self, method="shreve"):
        for node in self.nodes:
            if len(node) == 5:  # exactly one arc connected
                if (
                    node[1] == self.network[node[4]][3]
                    and node[2] == self.network[node[4]][4]
                ):
                    sink = self.network[node[4]][0]
                    if method == "shreve":
                        self.network[sink][7] = shreve(
                            sink, self.network[sink][5], self.network, self.nodes
                        )
                    else:
                        self.network[sink][7] = strahler(
                            sink, self.network[sink][5], self.network, self.nodes
                        )

    def load_attributes(self):
        for node in self.nodes:
            node_proj = merc_pro(*node[1:3], inverse=True)

            node_elevation = next(self.dem.sample([node_proj]))[0]
            node_runoff = next(self.flow.sample([node_proj]))[0]
            node_flow_acc = next(self.flowacc.sample([node_proj]))[0]

            node[3] = {
                "elevation": node_elevation,
                "runoff": node_runoff,
                "flow_acc": node_flow_acc,
            }

    def network_to_gdf(self):
        network_df = pd.DataFrame(
            columns=[
                "idx",
                "xs",
                "ys",
                "xe",
                "ye",
                "node_start",
                "node_end",
                "so",
                "length",
            ],
            data=self.network,
        )
        self.rivers_out = self.rivers.merge(
            network_df, how="left", left_index=True, right_index=True
        )
        self.rivers_out = self.rivers_out.to_crs(wgs84)

        nodes_for_df = [
            node[0:3] for node in self.nodes
        ]  # drop the extra columns that will confuse a df
        for node in self.nodes:
            nodes_for_df[node[0]].extend(list(node[3].values()))
            nodes_for_df[node[0]].append(node[4:])
        nodes_df = pd.DataFrame(
            columns=["idx", "x", "y", "elevation", "runoff", "flow_acc", "arcs"],
            data=nodes_for_df,
        )
        nodes_geometry = [Point(xy) for xy in zip(nodes_df["x"], nodes_df["y"])]
        self.nodes_out = gpd.GeoDataFrame(
            nodes_df, crs=self.rivers.crs, geometry=nodes_geometry
        )
        self.nodes_out = self.nodes_out.to_crs(wgs84)

    def local_flowacc(self):
        """
        For each node, subtract the flow accumulation for all upstream nodes, so that
        we calculate it's 'self discharge' only for it's directly contributing area.
        If anything goes negative, give its original flow_acc back.
        """

        self.nodes_out["flow_acc_local"] = 0
        for index, node in self.nodes_out.iterrows():
            actual_flow_acc = self.nodes_out.loc[index, "flow_acc"]
            for arc in node["arcs"]:
                if self.network[arc][6] == index:
                    subtract = self.nodes_out.loc[self.network[arc][5], "flow_acc"]
                    actual_flow_acc -= subtract
            self.nodes_out.loc[index, "flow_acc_local"] = actual_flow_acc
        self.nodes_out.loc[
            self.nodes_out["flow_acc_local"] < 0, "flow_acc_local"
        ] = self.nodes_out["flow_acc"]

    def calculate_runoff(self, kc_path):
        self.nodes_out, self.runoff_df = rainfall_runoff(self.nodes_out, kc_path)

    def calculate_dischange(self):
        self.nodes_out, self.rivers_out = discharge(
            self.nodes_out, self.network, self.rivers_out, self.flowacc, self.runoff_df,
        )

    def calculate_hydro(self, interval=1000, head_distance=500, get_range=False):
        """
        Start at the top of each river, get source GSCDxFlowAccxsquare size.
        Calculates flow rate in m3/s, head in m and power in Watts.

        Parameters
        ----------
        interval : int, optional (default 1000.)
            Interval in metres between testing points.
        head_distance : int, optional (default 500.)
            How far upstream to look to calculate head.
        get_range : boolean, optional (default False.)
            Whether to attempt to calculate mean, median, max etc.
        """

        # Choose the point_interval in metres
        # Loops through each stream and adds points down it at the specified interval
        # And creates a dict with these new point geometries

        # At 20 deg South:
        # 1 degree = 110704 metres -> 1 minute = 1845 metres -> 1 second = 30.75 metres
        # River network is 15 second resolution = 461 metres
        # Therefore each up_cell size is
        cell_area = (110704 / 60 / 60 * 15) ** 2
        rivers = self.rivers.to_crs(merc)

        rho = 998.57  # density of water, kg/m3
        g = 9.81  # acceleration due to gravity, m/2
        n = 0.5  # overall efficiency of system
        eta_t = 0.88
        eta_g = 0.96
        conv = 0.6

        hydro_points_dicts = []
        count = 0
        for _, row in rivers.iterrows():
            geom = shape(row["geometry"])
            length = geom.length
            for _, distance in enumerate(range(0, int(length), interval)):
                arcid = row["arcid"]
                up_cells = row["up_cells"]

                point = Point(
                    merc_pro(*list(geom.interpolate(distance).coords)[0], inverse=True)
                )
                upstream_point = Point(
                    merc_pro(
                        *list(geom.interpolate(distance - head_distance).coords)[0],
                        inverse=True
                    )
                )

                elevation = next(self.dem.sample(list(point.coords))).tolist()[0]
                upstream_elevation = next(
                    self.dem.sample(list(upstream_point.coords))
                ).tolist()[0]
                head = upstream_elevation - elevation

                runoff = next(self.flow.sample(list(point.coords))).tolist()[0]
                flowrate = (
                    runoff * up_cells * cell_area * (1 / 1000) * (1 / (8760 * 3600))
                )  # convert mm/year to m3/s

                power = rho * g * n * flowrate * head

                if head > 0 and flowrate > 0:
                    hydro_points_dicts.append(
                        {
                            "arcid": arcid,
                            "elevation": elevation,
                            "head": head,
                            "flowrate": flowrate,
                            "power": power,
                            "geometry": point,
                        }
                    )

                count += 1

        self.hydro_points = gpd.GeoDataFrame(hydro_points_dicts)
        self.hydro_points.crs = wgs84

        if get_range:
            # transfer discharge values from the river to the points
            self.hydro_points["discharge_accum"] = -99
            self.hydro_points["discharge_max"] = -99
            self.hydro_points["discharge_mean"] = -99
            self.hydro_points["discharge_min"] = -99
            for index, row in self.nodes_out.iterrows():
                self.hydro_points.loc[index, "discharge_accum"] = rivers.loc[
                    self.nodes_out.loc[index, "arcid"], "discharge_accum"
                ]
                self.hydro_points.loc[index, "discharge_max"] = rivers.loc[
                    self.nodes_out.loc[index, "arcid"], "discharge_max"
                ]
                self.hydro_points.loc[index, "discharge_mean"] = rivers.loc[
                    self.nodes_out.loc[index, "arcid"], "discharge_mean"
                ]
                self.hydro_points.loc[index, "discharge_min"] = rivers.loc[
                    self.nodes_out.loc[index, "arcid"], "discharge_min"
                ]

            pow_mult = rho * g * eta_t * eta_g * conv
            self.hydro_points["power_accum"] = (
                pow_mult * self.nodes_out["discharge_accum"] * self.nodes_out["head"]
            )
            self.hydro_points["power_max"] = (
                pow_mult * self.nodes_out["discharge_max"] * self.nodes_out["head"]
            )
            self.hydro_points["power_mean"] = (
                pow_mult * self.nodes_out["discharge_mean"] * self.nodes_out["head"]
            )
            self.hydro_points["power_min"] = (
                pow_mult * self.nodes_out["discharge_min"] * self.nodes_out["head"]
            )
