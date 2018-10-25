from __future__ import absolute_import, division, print_function, unicode_literals
import arcpy
import numpy as np
import pandas as pd
import heapq
from os.path import join
from datetime import datetime


field_flow_acc = 'flow_acc'
field_gscd = 'gscd'
field_flow_acc_local = 'flow_acc_local'
field_discharge_local = 'discharge_local'
field_discharge_accumulated = 'discharge'
field_discharge_max = 'dis_max'
field_discharge_mean = 'dis_mean'
field_discharge_min = 'dis_min'
field_arcs = 'arcs'
field_shreve = 'shreve'
field_land_type = 'land_type'
field_elevation = 'elevation'
field_head = 'head'
field_power = 'power'
field_power_max = 'power_max'
field_power_mean = 'power_mean'
field_power_min = 'power_min'
field_arc_id = 'arc_id'
field_et_ref = 'et_ref'
field_precip = 'precip_'
days_per_month = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30, 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}


def shreve(arc_index, direction_node_id, network, nodes):
    """
    Caclulate Shreve stream order instead of Strahler
    This ensures that a downstream river is always higher
    """
    up_stream_orders = []
    if len(nodes[direction_node_id]) == 4:
        network[arc_index][7] = 1
    else:
        for index, arc in enumerate(nodes[direction_node_id]):
            if index >= 3:
                if network[arc][0] != arc_index:
                    if network[arc][5] != direction_node_id:
                        up_stream_orders.append(shreve(arc, network[arc][5], network, nodes))
                    else:
                        up_stream_orders.append(shreve(arc, network[arc][6], network, nodes))

        max_orders = heapq.nlargest(2, up_stream_orders)
        if len(max_orders) == 2:
            order = 0 + max_orders[0] + max_orders[1]
        else:
            order = 0 + max(up_stream_orders)

        network[arc_index][7] = order

    print('so {}'.format(arc_index))
    return network[arc_index][7]


def strahler(arc_index, direction_node_id, network, nodes):
    """
    This function is nearly verbatim from Gleyzer2004 algorithm
    But excludes the code to create river segments

    :param arc_index: the index of an arc in the networks array
    :param direction_node_id: the index of that arc's upstream node
    :param network: the network of arcs
    :param nodes: contains all the nodes and their connections
    """
    up_stream_orders = []
    if len(nodes[direction_node_id]) == 1:
        network[arc_index][7] = 1
    else:
        for arc in nodes[direction_node_id]:
            if network[arc][0] != arc_index:
                if network[arc][5] != direction_node_id:
                    up_stream_orders.append(strahler(arc, network[arc][5], network, nodes))
                else:
                    up_stream_orders.append(strahler(arc, network[arc][6], network, nodes))
        max_order = 0
        max_order_count = 0
        for order in up_stream_orders:
            if order > max_order:
                max_order = order
                max_order_count = 1
            elif order == max_order:
                max_order_count += 1
        if max_order_count > 1:
            network[arc_index][7] = max_order + 1
        else:
            network[arc_index][7] = max_order
    print('so {}'.format(arc_index))
    return network[arc_index][7]


class HydroModeller:
    def __init__(self, fc_rivers, fc_flow_acc, fc_gscd, fc_land_type, fc_elevation, fc_points, fc_et_ref, fc_nodes,
                 fc_precip_prefix):
        self.fc_rivers = fc_rivers
        self.fc_flow_acc = fc_flow_acc
        self.fc_gscd = fc_gscd
        self.fc_land_type = fc_land_type
        self.fc_elevation = fc_elevation
        self.fc_points = fc_points
        self.fc_et_ref = fc_et_ref
        self.fc_nodes = fc_nodes
        self.fc_precip_prefix = fc_precip_prefix

        self.network = None
        self.nodes = None
        self.network_df = None
        self.nodes_df = None
        # self.precipitation_df = None
        self.runoff_df = None
        self.discharge_df = None
        self.points_df = None
        self.precipitation_df = None

    def create_network(self):
        """
        Create the network and nodes from the HydroSHEDS layer, ready to be stream ordered

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
        3..   arc indices...
        """

        # import the start and end end points into the network array of arcs
        # and simultaneously import them into a separate node array (so two nodes per arc)
        # there will be duplicate nodes, but they duplicates are ignored later
        counter = 0
        length = int(arcpy.GetCount_management(self.fc_rivers).getOutput(0))
        self.network = np.empty([length, 9], dtype=np.int32)
        self.nodes = []
        with arcpy.da.SearchCursor(self.fc_rivers, ['SHAPE@']) as cursor:
            for row in cursor:
                # the geometry variables (presumably int so that they actual match when comparing)
                fx = int(row[0].firstPoint.X)
                fy = int(row[0].firstPoint.Y)
                lx = int(row[0].lastPoint.X)
                ly = int(row[0].lastPoint.Y)

                # add arc_length as a determinant of how much water contributes downstream
                # not sure if it's a good idea going forward...
                arc_length = int(row[0].length)

                # store the index as a column for easier access
                # the last column is for stream order
                self.network[counter] = [counter, fx, fy, lx, ly, -99, -99, -99, arc_length]

                # create nodes for the start point and end point separately
                self.nodes.append([2*counter, fx, fy])
                self.nodes.append([2*counter+1, lx, ly])

                counter += 1
                print('imported {}'.format(counter))

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
            print('indexed {}'.format(arc[0]))

        # for every node, add references to every arc that connects to it
        for arc in self.network:
            # tell the arc's starting node that it exists
            self.nodes[arc[5]].append(arc[0])
            # and tell the arc's ending node that it exists
            self.nodes[arc[6]].append(arc[0])
            print('nodes {}'.format(arc[0]))

    def calc_stream_order(self):
        """
        """

        for node in self.nodes:
            if len(node) == 4:  # only one arc connected
                if node[1] == self.network[node[3]][3] and node[2] == self.network[node[3]][4]:
                    sink = self.network[node[3]][0]
                    self.network[sink][7] = shreve(sink, self.network[sink][5], self.network, self.nodes)

    def prepare_nodes(self):
        """
        Get the geospatial characteristics at each node
        """

        # Create a points feature class, and add points using the coordinates of every node
        desc = arcpy.Describe(self.fc_rivers)
        nodes_arr = np.empty(len(self.nodes), dtype=[(str('x'), np.float32), (str('y'), np.float32)])
        nodes_arr['x'] = [n[1] for n in self.nodes]
        nodes_arr['y'] = [n[2] for n in self.nodes]

        if arcpy.Exists(join(arcpy.env.workspace, self.fc_nodes)):
            arcpy.Delete_management(join(arcpy.env.workspace, self.fc_nodes))
        arcpy.da.NumPyArrayToFeatureClass(nodes_arr, join(arcpy.env.workspace, self.fc_nodes), ('x', 'y'),
                                          spatial_reference=desc.spatialReference)

        # Extract the Flow accumulation and gscd values into these points, from the rasters
        arcpy.CheckOutExtension('Spatial')
        arcpy.sa.ExtractMultiValuesToPoints(self.fc_nodes, [[self.fc_flow_acc, field_flow_acc]])
        arcpy.sa.ExtractMultiValuesToPoints(self.fc_nodes, [[self.fc_gscd, field_gscd]])
        arcpy.sa.ExtractMultiValuesToPoints(self.fc_nodes, [[self.fc_land_type, field_land_type]])
        arcpy.sa.ExtractMultiValuesToPoints(self.fc_nodes, [[self.fc_et_ref, field_et_ref]])

        for i in range(1, 13):
            arcpy.sa.ExtractMultiValuesToPoints(self.fc_nodes, [['{}{}'.format(self.fc_precip_prefix, i),
                                                                 '{}{}'.format(field_precip, i)]])
        arcpy.CheckInExtension('Spatial')

        # load flowacc and runoff from the shapefile we just created
        self.nodes_df = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(self.fc_nodes, [
            'SHAPE@X', 'SHAPE@Y', field_flow_acc, field_gscd, field_land_type,
            field_et_ref] + ['{}{}'.format(field_precip, i) for i in range(1,13)], skip_nulls=False, null_value=0))
        arc_refs = pd.Series(dtype=tuple, index=self.nodes_df.index.tolist())
        for index, node in self.nodes_df.iterrows():
            arc_refs[index] = self.nodes[index][3:]
        self.nodes_df[field_arcs] = arc_refs

    def calc_flow_acc(self):
        """
        """

        # For each node, subtract the flow accumulation for all upstream nodes,
        # so that we calculate it's 'self discharge' only for it's directly contributing area
        # if anything goes negative, give its original flow_acc back
        self.nodes_df[field_flow_acc_local] = 0
        for index, node in self.nodes_df.iterrows():
            actual_flowacc = self.nodes_df[field_flow_acc][index]
            for arc in node[field_arcs]:
                if self.network[arc][6] == index:  # the up arc flows *into* the node
                    subtract = self.nodes_df[field_flow_acc][self.network[arc][5]]
                    actual_flowacc -= subtract  # subtract upstream flowAcc
            self.nodes_df.loc[index, field_flow_acc_local] = actual_flowacc
        self.nodes_df.loc[self.nodes_df[field_flow_acc_local] < 0, field_flow_acc_local] = self.nodes_df[field_flow_acc]

    def rainfall_runoff(self, default_precip_effectiveness, default_runoff_to_gw_fraction, runoff_calibration_accuracy):
        """
        Calculate the runoff from the rainfall using the rainfall-runoff method
        Calibrate against gscd annual values
        """
        # precipitation_se = pd.read_excel('precip.xlsx', index_col=0, squeeze=True)
        # self.precipitation_df = pd.DataFrame(columns=precipitation_se.index, index=self.nodes_df.index)
        # for index, row in self.precipitation_df.iterrows():
        #    self.precipitation_df.loc[index] = precipitation_se

        self.runoff_df = pd.DataFrame(columns=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], index=self.nodes_df.index)

        k_c = pd.read_excel('k_c.xlsx', index_col=0)

        # formula from http://www.weap21.org/webhelp/hydrology.htm
        # for each point, calculate the runoff for each month, then compare to gscd annual value and modify params
        self.nodes_df['precip_effective'] = default_precip_effectiveness  # default starting values
        self.nodes_df['runoff_to_gw_fraction'] = default_runoff_to_gw_fraction
        for index, row in self.runoff_df.iterrows():
            print('calibrate {}'.format(index))
            counter = 0
            while True:
                for col in self.runoff_df:
                    et_ref = self.nodes_df.loc[index, field_et_ref]
                    # precip = self.precipitation_df.loc[index, col]
                    precip = self.nodes_df.loc[index, '{}{}'.format(field_precip, col)] * \
                             3600 * 24 * days_per_month[col]
                    precip_avail_for_et = precip * self.nodes_df.loc[index, 'precip_effective']
                    et_potential = et_ref * k_c.loc[col, self.nodes_df.loc[index, field_land_type]]
                    runoff = max(0, precip_avail_for_et-et_potential) + \
                        precip * (1-self.nodes_df.loc[index, 'precip_effective'])
                    runoff_to_surface = runoff * (1 - self.nodes_df.loc[index, 'runoff_to_gw_fraction'])

                    self.runoff_df.loc[index, col] = runoff_to_surface

                ref_value = self.nodes_df.loc[index, field_gscd]
                calc_value = self.runoff_df.loc[index].mean() * 12  # to compare with annual gscd

                if counter > 10:
                    break
                elif abs((calc_value - ref_value) / ref_value) < runoff_calibration_accuracy:
                    break
                else:
                    counter += 1
                    self.nodes_df.loc[index, 'precip_effective'] *= 1 + sorted([-0.5, (calc_value - ref_value) / ref_value / 10, 0.5])[1]
                    self.nodes_df.loc[index, 'runoff_to_gw_fraction'] *= 1 + sorted([-0.5, (calc_value - ref_value) / ref_value / 10, 0.5])[1]

    def calc_discharge(self, water_loss_factor, max_water_loss_fraction):
        """
        """

        area = float(arcpy.GetRasterProperties_management(self.fc_flow_acc, 'CELLSIZEX').getOutput(0)) ** 2

        # calculate self discharge for each node
        # Q [m3/s] = runoff [mm] * flowacc [number] * area [m2] / 8760 [h/year] * 3600 [s/h] * 1000 [mm/m]
        self.discharge_df = pd.DataFrame(columns=self.runoff_df.columns, index=self.runoff_df.index)
        for index, row in self.discharge_df.iterrows():
            print('discharge {}'.format(index))
            for col in self.discharge_df:
                self.discharge_df.loc[index, col] = \
                    self.runoff_df.loc[index, col]*self.nodes_df.loc[index, field_flow_acc_local] * \
                    area / (24*days_per_month[col]*3600*1000)

        self.nodes_df[field_discharge_local] = \
            self.nodes_df[field_gscd] * self.nodes_df[field_flow_acc_local] * area / (8760*3600*1000)
        self.nodes_df[field_discharge_accumulated] = self.nodes_df[field_discharge_local]

        # start from Shreve order 1, and contribute all discharges from upstream nodes to downstream nodes
        for so in range(1, max(self.network.T[7])+1):
            for arc in self.network:
                if arc[7] == so:
                    print('contribute {}'.format(arc[0]))
                    self.nodes_df.loc[arc[6], field_discharge_accumulated] += \
                        self.nodes_df.loc[arc[5], field_discharge_accumulated] * \
                        (1 - max(max_water_loss_fraction, water_loss_factor * arc[8]))
                    for col in self.discharge_df:
                        self.discharge_df.loc[arc[6], col] += \
                            self.discharge_df.loc[arc[5], col] * (1 - max(max_water_loss_fraction,
                                                                          water_loss_factor * arc[8]))

        self.nodes_df[field_discharge_max] = -99
        self.nodes_df[field_discharge_mean] = -99
        self.nodes_df[field_discharge_min] = -99
        for index, node in self.nodes_df.iterrows():
            self.nodes_df.loc[index, field_discharge_max] = self.discharge_df.loc[index].max()
            self.nodes_df.loc[index, field_discharge_mean] = self.discharge_df.loc[index].mean()
            self.nodes_df.loc[index, field_discharge_min] = self.discharge_df.loc[index].min()

        # add the discharge from the upstream node of each arc to that arc
        self.network_df = pd.DataFrame(self.network, columns=[field_arc_id, 'fx', 'fy', 'lx', 'ly',
                                                              'nif', 'nil', field_shreve, 'arc_length'])
        self.network_df[field_discharge_accumulated] = -99
        self.network_df[field_discharge_max] = -99
        self.network_df[field_discharge_mean] = -99
        self.network_df[field_discharge_min] = -99
        for index, arc in self.network_df.iterrows():
            self.network_df.loc[index, field_discharge_accumulated] = \
                self.nodes_df[field_discharge_accumulated][arc['nif']]
            self.network_df.loc[index, field_discharge_max] = self.nodes_df[field_discharge_max][arc['nif']]
            self.network_df.loc[index, field_discharge_mean] = self.nodes_df[field_discharge_mean][arc['nif']]
            self.network_df.loc[index, field_discharge_min] = self.nodes_df[field_discharge_min][arc['nif']]

        # add the results back into the ArcGIS feature
        arcpy.AddField_management(self.fc_rivers, field_arc_id, 'LONG')
        arcpy.AddField_management(self.fc_rivers, field_shreve, 'LONG')
        arcpy.AddField_management(self.fc_rivers, field_discharge_accumulated, 'DOUBLE')
        arcpy.AddField_management(self.fc_rivers, field_discharge_max, 'DOUBLE')
        arcpy.AddField_management(self.fc_rivers, field_discharge_mean, 'DOUBLE')
        arcpy.AddField_management(self.fc_rivers, field_discharge_min, 'DOUBLE')
        with arcpy.da.UpdateCursor(self.fc_rivers, [field_arc_id, field_shreve, field_discharge_accumulated,
                                                    field_discharge_max, field_discharge_mean,
                                                    field_discharge_min]) as cursor:
            for index, row in enumerate(cursor, start=0):
                row[0] = self.network_df.loc[index, field_arc_id]
                row[1] = self.network_df.loc[index, field_shreve]
                row[2] = self.network_df.loc[index, field_discharge_accumulated]
                row[3] = self.network_df.loc[index, field_discharge_max]
                row[4] = self.network_df.loc[index, field_discharge_mean]
                row[5] = self.network_df.loc[index, field_discharge_min]
                cursor.updateRow(row)
                print('added {}'.format(index))

    def calc_hydro_potential(self, interval, eta_t, eta_g, conv):
        """
        Insert points at intervals along every stream and estimate the hydro potential at each
        """

        desc = arcpy.Describe(self.fc_rivers)
        arcpy.CreateFeatureclass_management(arcpy.env.workspace, self.fc_points, geometry_type='POINT',
                                            spatial_reference=desc.spatialReference)
        arcpy.AddField_management(self.fc_points, field_arc_id, 'LONG')

        with arcpy.da.SearchCursor(self.fc_rivers, ['SHAPE@', field_arc_id]) as search_cursor:
            with arcpy.da.InsertCursor(self.fc_points, ['SHAPE@', field_arc_id]) as insert_cursor:
                for row in search_cursor:
                    print('pointsifying {}'.format(row[1]))
                    line = row[0]

                    if line:
                        cur_length = interval
                        max_position = line.length
                        insert_cursor.insertRow([line.firstPoint, row[1]])
                        while cur_length < max_position:
                            insert_cursor.insertRow([line.positionAlongLine(cur_length, False), row[1]])
                            cur_length += interval

        arcpy.CheckOutExtension('Spatial')
        arcpy.sa.ExtractMultiValuesToPoints(self.fc_points, [[self.fc_elevation, field_elevation]])
        arcpy.CheckInExtension('Spatial')

        self.points_df = pd.DataFrame(arcpy.da.TableToNumPyArray(self.fc_points, (field_arc_id, field_elevation),
                                                            skip_nulls=False, null_value=0))

        # transfer discharge values from the river to the points
        self.points_df[field_discharge_accumulated] = -99
        self.points_df[field_discharge_max] = -99
        self.points_df[field_discharge_mean] = -99
        self.points_df[field_discharge_min] = -99
        for index, row in self.points_df.iterrows():
            self.points_df.loc[index, field_discharge_accumulated] = \
                self.network_df.loc[self.points_df.loc[index, field_arc_id], field_discharge_accumulated]
            self.points_df.loc[index, field_discharge_max] = \
                self.network_df.loc[self.points_df.loc[index, field_arc_id], field_discharge_max]
            self.points_df.loc[index, field_discharge_mean] = \
                self.network_df.loc[self.points_df.loc[index, field_arc_id], field_discharge_mean]
            self.points_df.loc[index, field_discharge_min] = \
                self.network_df.loc[self.points_df.loc[index, field_arc_id], field_discharge_min]

        # calculate the elevation difference (head) between successive points
        self.points_df[field_head] = -99
        prev_fid = 0
        prev_elev = 0
        for index, point in self.points_df.iterrows():
            print('head {}'.format(index))
            if point[field_elevation]:  # check to make sure an elevation entry exists
                current_fid = point[field_arc_id]
                current_elev = point[field_elevation]

                if (current_fid == prev_fid) and ((prev_elev - current_elev) > 0):
                    self.points_df.loc[index, field_head] = prev_elev - current_elev
                else:
                    self.points_df.loc[index, field_head] = 0

                prev_fid = current_fid
                prev_elev = current_elev
            else:
                self.points_df.loc[index, field_head] = 0
                prev_fid = point[field_arc_id]
                prev_elev = 0

        # calculate the power in watts based on Alex's formula
        # P = rho * g * eta_t * eta_g * conv * Q * deltaH
        rho = 1000  # density
        g = 9.80665  # gravity
        self.points_df[field_power] = rho * g * eta_t * eta_g * conv * \
            self.points_df[field_discharge_accumulated] * self.points_df[field_head]
        self.points_df[field_power_max] = rho * g * eta_t * eta_g * conv * \
            self.points_df[field_discharge_max] * self.points_df[field_head]
        self.points_df[field_power_mean] = rho * g * eta_t * eta_g * conv * \
            self.points_df[field_discharge_mean] * self.points_df[field_head]
        self.points_df[field_power_min] = rho * g * eta_t * eta_g * conv * \
            self.points_df[field_discharge_min] * self.points_df[field_head]

        # add the results back into the ArcGIS feature
        arcpy.AddField_management(self.fc_points, field_discharge_accumulated, 'DOUBLE')
        arcpy.AddField_management(self.fc_points, field_discharge_max, 'DOUBLE')
        arcpy.AddField_management(self.fc_points, field_discharge_mean, 'DOUBLE')
        arcpy.AddField_management(self.fc_points, field_discharge_min, 'DOUBLE')
        arcpy.AddField_management(self.fc_points, field_head, 'FLOAT')
        arcpy.AddField_management(self.fc_points, field_power, 'DOUBLE')
        arcpy.AddField_management(self.fc_points, field_power_max, 'DOUBLE')
        arcpy.AddField_management(self.fc_points, field_power_mean, 'DOUBLE')
        arcpy.AddField_management(self.fc_points, field_power_min, 'DOUBLE')
        with arcpy.da.UpdateCursor(self.fc_points, [field_discharge_accumulated, field_discharge_max,
                                                    field_discharge_mean, field_discharge_min, field_head, field_power,
                                                    field_power_max, field_power_mean, field_power_min]) as cursor:
            for index, row in enumerate(cursor, start=0):
                row[0] = self.points_df.loc[index, field_discharge_accumulated]
                row[1] = self.points_df.loc[index, field_discharge_max]
                row[2] = self.points_df.loc[index, field_discharge_mean]
                row[3] = self.points_df.loc[index, field_discharge_min]
                row[4] = self.points_df.loc[index, field_head]
                row[5] = self.points_df.loc[index, field_power]
                row[6] = self.points_df.loc[index, field_power_max]
                row[7] = self.points_df.loc[index, field_power_mean]
                row[8] = self.points_df.loc[index, field_power_min]
                cursor.updateRow(row)
                print('added {}'.format(index))


def runner():
    arcpy.env.workspace = r'C:/Users/Chris/Documents/GIS/HydroModeller.gdb'
    arcpy.env.overwriteOutput = True
    arcpy.env.addOutputsToMap = False

    # TODO maybe need to aggregate elevation, as doesn't capture always
    modeller = HydroModeller(fc_rivers='rivers_ug_small_projected',
                             fc_flow_acc='flow_acc_af_projected',
                             fc_gscd='gscd_qmean_projected',
                             fc_land_type='land_type_projected',
                             fc_elevation='elevation_af_projected',
                             fc_points='hydro_points_ug2',
                             fc_et_ref='et_ref_projected',
                             fc_nodes='_temp_nodes2',
                             fc_precip_prefix='ave_p_')

    modeller.create_network()
    modeller.calc_stream_order()
    modeller.prepare_nodes()
    modeller.calc_flow_acc()
    modeller.rainfall_runoff(default_precip_effectiveness=0.5,
                             default_runoff_to_gw_fraction=0.5,
                             runoff_calibration_accuracy=0.2)
    modeller.calc_discharge(water_loss_factor=0.00001, max_water_loss_fraction=0.2)
    modeller.calc_hydro_potential(interval=1000, eta_t=0.88, eta_g=0.96, conv=0.6)


if __name__ == "__main__":
    startTime = datetime.now()
    runner()
    print('time {}'.format(datetime.now() - startTime))
