from __future__ import absolute_import, division, print_function, unicode_literals
from os import path
arcpy = path.join(r'C:\Program Files (x86)\ArcGIS\Desktop10.5\arcpy\arcpy', "arcpy")
import arcpy
import numpy as np
import pandas as pd
import heapq
from os.path import join
from datetime import datetime

nodes_fc = '_temp_nodes'
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


def create_network(fc_rivers):
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
    length = int(arcpy.GetCount_management(fc_rivers).getOutput(0))
    network = np.empty([length, 8], dtype=np.int32)
    nodes = []
    with arcpy.da.SearchCursor(fc_rivers, ['SHAPE@']) as cursor:
        for row in cursor:
            # the geometry variables (presumably int so that they actual match when comparing)
            fx = int(row[0].firstPoint.X)
            fy = int(row[0].firstPoint.Y)
            lx = int(row[0].lastPoint.X)
            ly = int(row[0].lastPoint.Y)

            # store the index as a column for easier access
            # the last column is for stream order
            network[counter] = [counter, fx, fy, lx, ly, -99, -99, -99]

            # create nodes for the start point and end point separately
            nodes.append([2*counter, fx, fy])
            nodes.append([2*counter+1, lx, ly])

            counter += 1
            print('imported {}'.format(counter))

    # the most time consuming part of the script
    # runs through every arc and every node and where there is a match, add
    #  a reference to the index into the network array
    # once both nodes have been found for an arc, break that for loop and
    #  start with the next arc
    for arc in network:
        match_f, match_l = False, False
        for node in nodes:
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
    for arc in network:
        # tell the arc's starting node that it exists
        nodes[arc[5]].append(arc[0])
        # and tell the arc's ending node that it exists
        nodes[arc[6]].append(arc[0])
        print('nodes {}'.format(arc[0]))

    return network, nodes


def calc_stream_order(network, nodes):
    """
    """

    for node in nodes:
        if len(node) == 4:  # only one arc connected
            if node[1] == network[node[3]][3] and node[2] == network[node[3]][4]:
                sink = network[node[3]][0]
                network[sink][7] = shreve(sink, network[sink][5], network, nodes)

    return network


def prepare_nodes(fc_rivers, nodes, fc_flow_acc, fc_gscd, fc_land_type, fc_et_ref):
    """
    Get the geospatial characteristics at each node
    """

    # Create a points feature class, and add points using the coordinates of every node
    desc = arcpy.Describe(fc_rivers)
    nodes_arr = np.empty(len(nodes), dtype=[(str('x'), np.float32), (str('y'), np.float32)])
    nodes_arr['x'] = [n[1] for n in nodes]
    nodes_arr['y'] = [n[2] for n in nodes]

    if arcpy.Exists(join(arcpy.env.workspace, nodes_fc)):
        arcpy.Delete_management(join(arcpy.env.workspace, nodes_fc))
    arcpy.da.NumPyArrayToFeatureClass(nodes_arr, join(arcpy.env.workspace, nodes_fc), ('x', 'y'),
                                      spatial_reference=desc.spatialReference)

    # Extract the Flow accumulation and gscd values into these points, from the rasters
    arcpy.CheckOutExtension('Spatial')
    arcpy.sa.ExtractMultiValuesToPoints(nodes_fc, [[fc_flow_acc, field_flow_acc]])
    arcpy.sa.ExtractMultiValuesToPoints(nodes_fc, [[fc_gscd, field_gscd]])
    arcpy.sa.ExtractMultiValuesToPoints(nodes_fc, [[fc_land_type, field_land_type]])
    arcpy.sa.ExtractMultiValuesToPoints(nodes_fc, [[fc_et_ref, field_et_ref]])

    # precipitation should be extracted to a separate df?
    # arcpy.sa.ExtractMultiValuesToPoints(nodes_fc, [[fc_precipitation, field_precipitation]])
    arcpy.CheckInExtension('Spatial')

    # load flowacc and runoff from the shapefile we just created
    nodes_df = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(nodes_fc, ('SHAPE@X', 'SHAPE@Y', field_flow_acc,
                                                                         field_gscd, field_land_type, field_et_ref),
                                                              skip_nulls=False, null_value=0))
    arc_refs = pd.Series(dtype=tuple, index=nodes_df.index.tolist())
    for index, node in nodes_df.iterrows():
        arc_refs[index] = nodes[index][3:]
    nodes_df[field_arcs] = arc_refs

    return nodes_df


def calc_flow_acc(nodes_df, network):
    """
    """

    # For each node, subtract the flow accumulation for all upstream nodes,
    # so that we calculate it's 'self discharge' only for it's directly contributing area
    # if anything goes negative, give its original flow_acc back
    nodes_df[field_flow_acc_local] = 0
    for index, node in nodes_df.iterrows():
        actual_flowacc = nodes_df[field_flow_acc][index]
        for arc in node[field_arcs]:
            if network[arc][6] == index:  # the up arc flows *into* the node
                subtract = nodes_df[field_flow_acc][network[arc][5]]
                actual_flowacc -= subtract  # subtract upstream flowAcc
        nodes_df.loc[index, field_flow_acc_local] = actual_flowacc
    nodes_df.loc[nodes_df[field_flow_acc_local] < 0, field_flow_acc_local] = nodes_df[field_flow_acc]

    return nodes_df


def rainfall_runoff(nodes_df):
    """
    Calculate the runoff from the rainfall using the rainfall-runoff method
    Calibrate against gscd annual values
    """
    precipitation_se = pd.read_excel('precip.xlsx', index_col=0, squeeze=True)
    precipitation_df = pd.DataFrame(columns=precipitation_se.index, index=nodes_df.index)
    for index, row in precipitation_df.iterrows():
        precipitation_df.loc[index] = precipitation_se

    runoff_df = pd.DataFrame(columns=precipitation_df.columns, index=precipitation_df.index)

    # k_c should vary by month...
    k_c = {0: 0.6,  # there shouldn't be a zero, this is error handling
           1: 0.6,
           2: 0.9,
           3: 0.6,
           4: 0.6,
           5: 0.6,
           6: 0.5,
           7: 0.6,
           8: 0.6,
           9: 0.6,
           10: 0.1,
           11: 0.6}

    # formula from http://www.weap21.org/webhelp/hydrology.htm
    # for each point, calculate the runoff for each month, then compare to gscd annual value and modify params if need
    nodes_df['precip_effective'] = 0.5  # default starting values
    nodes_df['runoff_to_gw_fraction'] = 0.5
    for index, row in runoff_df.iterrows():
        print('calibrate {}'.format(index))
        counter = 0
        while True:
            for col in runoff_df:
                et_ref = nodes_df.loc[index, field_et_ref]
                precip = precipitation_df.loc[index, col]
                precip_avail_for_et = precip * nodes_df.loc[index, 'precip_effective']
                et_potential = et_ref * k_c[nodes_df.loc[index, field_land_type]]
                runoff = max(0, precip_avail_for_et-et_potential) + precip * (1-nodes_df.loc[index, 'precip_effective'])
                runoff_to_surface = runoff * (1 - nodes_df.loc[index, 'runoff_to_gw_fraction'])

                runoff_df.loc[index, col] = runoff_to_surface

            ref_value = nodes_df.loc[index, field_gscd]
            calc_value = runoff_df.loc[index].mean() * 12  # to compare with annual gscd

            if counter > 10:
                break
            elif abs((calc_value - ref_value) / ref_value) < 0.5:
                break
            else:
                counter += 1
                nodes_df.loc[index, 'precip_effective'] *= 1 + (calc_value - ref_value) / ref_value / 10
                nodes_df.loc[index, 'runoff_to_gw_fraction'] *= 1 + (calc_value - ref_value) / ref_value / 10

    return nodes_df, runoff_df


def calc_discharge(nodes_df, runoff_df, network, fc_flow_acc, fc_rivers):
    """
    """

    area = float(arcpy.GetRasterProperties_management(fc_flow_acc, 'CELLSIZEX').getOutput(0)) ** 2

    # calculate self discharge for each node
    # Q [m3/s] = runoff [mm] * flowacc [number] * area [m2] / 8760 [h/year] * 3600 [s/h] * 1000 [mm/m]
    discharge_df = pd.DataFrame(columns=runoff_df.columns, index=runoff_df.index)
    for index, row in discharge_df.iterrows():
        print('discharge {}'.format(index))
        for col in discharge_df:
            # TODO need to take into account different month lengths
            discharge_df.loc[index, col] = runoff_df.loc[index, col]*nodes_df.loc[index, field_flow_acc_local]*area / (
                730*3600*1000)

    nodes_df[field_discharge_local] = nodes_df[field_gscd] * nodes_df[field_flow_acc_local] * area / (
        8760*3600*1000)
    nodes_df[field_discharge_accumulated] = nodes_df[field_discharge_local]

    # start from Shreve order 1, and contribute all discharges from upstream nodes to downstream nodes
    water_loss = 0.2
    for so in range(1, max(network.T[7])+1):
        for arc in network:
            if arc[7] == so:
                print('contribute {}'.format(arc[0]))
                nodes_df.loc[arc[6], field_discharge_accumulated] += nodes_df.loc[arc[5], field_discharge_accumulated]*(
                    1-water_loss)
                for col in discharge_df:
                    discharge_df.loc[arc[6], col] += discharge_df.loc[arc[5], col] * (1 - water_loss)

    nodes_df[field_discharge_max] = -99
    nodes_df[field_discharge_mean] = -99
    nodes_df[field_discharge_min] = -99
    for index, node in nodes_df.iterrows():
        nodes_df.loc[index, field_discharge_max] = discharge_df.loc[index].max()
        nodes_df.loc[index, field_discharge_mean] = discharge_df.loc[index].mean()
        nodes_df.loc[index, field_discharge_min] = discharge_df.loc[index].min()

    # add the discharge from the upstream node of each arc to that arc
    network_df = pd.DataFrame(network, columns=[field_arc_id, 'fx', 'fy', 'lx', 'ly', 'nif', 'nil', field_shreve])
    network_df[field_discharge_accumulated] = -99
    network_df[field_discharge_max] = -99
    network_df[field_discharge_mean] = -99
    network_df[field_discharge_min] = -99
    for index, arc in network_df.iterrows():
        network_df.loc[index, field_discharge_accumulated] = nodes_df[field_discharge_accumulated][arc['nif']]
        network_df.loc[index, field_discharge_max] = nodes_df[field_discharge_max][arc['nif']]
        network_df.loc[index, field_discharge_mean] = nodes_df[field_discharge_mean][arc['nif']]
        network_df.loc[index, field_discharge_min] = nodes_df[field_discharge_min][arc['nif']]

    # add the results back into the ArcGIS feature
    arcpy.AddField_management(fc_rivers, field_arc_id, 'LONG')
    arcpy.AddField_management(fc_rivers, field_shreve, 'LONG')
    arcpy.AddField_management(fc_rivers, field_discharge_accumulated, 'DOUBLE')
    arcpy.AddField_management(fc_rivers, field_discharge_max, 'DOUBLE')
    arcpy.AddField_management(fc_rivers, field_discharge_mean, 'DOUBLE')
    arcpy.AddField_management(fc_rivers, field_discharge_min, 'DOUBLE')
    with arcpy.da.UpdateCursor(fc_rivers, [field_arc_id, field_shreve, field_discharge_accumulated, field_discharge_max,
                                           field_discharge_mean, field_discharge_min]) as cursor:
        for index, row in enumerate(cursor, start=0):
            row[0] = network_df.loc[index, field_arc_id]
            row[1] = network_df.loc[index, field_shreve]
            row[2] = network_df.loc[index, field_discharge_accumulated]
            row[3] = network_df.loc[index, field_discharge_max]
            row[4] = network_df.loc[index, field_discharge_mean]
            row[5] = network_df.loc[index, field_discharge_min]
            cursor.updateRow(row)
            print('added {}'.format(index))

    return network_df, nodes_df, discharge_df


def calc_hydro_potential(fc_rivers, fc_elevation, interval, fc_points, network_df):
    """
    Insert points at intervals along every stream and estimate the hydro potential at each
    """

    desc = arcpy.Describe(fc_rivers)
    arcpy.CreateFeatureclass_management(arcpy.env.workspace, fc_points, geometry_type='POINT',
                                        spatial_reference=desc.spatialReference)
    arcpy.AddField_management(fc_points, field_arc_id, 'LONG')

    with arcpy.da.SearchCursor(fc_rivers, ['SHAPE@', field_arc_id]) as search_cursor:
        with arcpy.da.InsertCursor(fc_points, ['SHAPE@', field_arc_id]) as insert_cursor:
            for row in search_cursor:
                line = row[0]

                if line:
                    cur_length = interval
                    max_position = line.length
                    insert_cursor.insertRow([line.firstPoint, row[1]])
                    while cur_length < max_position:
                        insert_cursor.insertRow([line.positionAlongLine(cur_length, False), row[1]])
                        cur_length += interval

    arcpy.CheckOutExtension('Spatial')
    arcpy.sa.ExtractMultiValuesToPoints(fc_points, [[fc_elevation, field_elevation]])
    arcpy.CheckInExtension('Spatial')

    points_df = pd.DataFrame(arcpy.da.TableToNumPyArray(fc_points, (field_arc_id, field_elevation),
                                                        skip_nulls=False, null_value=0))

    points_df[field_discharge_accumulated] = -99
    points_df[field_discharge_max] = -99
    points_df[field_discharge_mean] = -99
    points_df[field_discharge_min] = -99
    for index, row in points_df.iterrows():
        points_df.loc[index, field_discharge_accumulated] = network_df.loc[points_df.loc[index, field_arc_id],
                                                                           field_discharge_accumulated]
        points_df.loc[index, field_discharge_max] = network_df.loc[points_df.loc[index, field_arc_id],
                                                                   field_discharge_max]
        points_df.loc[index, field_discharge_mean] = network_df.loc[points_df.loc[index, field_arc_id],
                                                                    field_discharge_mean]
        points_df.loc[index, field_discharge_min] = network_df.loc[points_df.loc[index, field_arc_id],
                                                                   field_discharge_min]

    # calculate the elevation difference (head) between successive points
    points_df[field_head] = -99
    prev_fid = 0
    prev_elev = 0
    for index, point in points_df.iterrows():
        if point[field_elevation]:  # check to make sure an elevation entry exists
            current_fid = point[field_arc_id]
            current_elev = point[field_elevation]

            if (current_fid == prev_fid) and ((prev_elev - current_elev) > 0):
                points_df.loc[index, field_head] = prev_elev - current_elev
            else:
                points_df.loc[index, field_head] = 0

            prev_fid = current_fid
            prev_elev = current_elev
        else:
            points_df.loc[index, field_head] = 0
            prev_fid = point[field_arc_id]
            prev_elev = 0

    # calculate the power in watts based on Alex's formula
    # P = rho * g * nt * ng * conv * Q * deltaH
    rho = 1000  # density
    g = 9.81  # gravity
    nt = 0.88  # turbine efficiency
    ng = 0.96  # generator efficiency
    conv = 0.6  # conversion factor for environmental flow deduction
    points_df[field_power] = rho * g * nt * ng * conv * points_df[field_discharge_accumulated] * points_df[field_head]
    points_df[field_power_max] = rho * g * nt * ng * conv * points_df[field_discharge_max] * points_df[field_head]
    points_df[field_power_mean] = rho * g * nt * ng * conv * points_df[field_discharge_mean] * points_df[field_head]
    points_df[field_power_min] = rho * g * nt * ng * conv * points_df[field_discharge_min] * points_df[field_head]

    # add the results back into the ArcGIS feature
    arcpy.AddField_management(fc_points, field_discharge_accumulated, 'DOUBLE')
    arcpy.AddField_management(fc_points, field_discharge_max, 'DOUBLE')
    arcpy.AddField_management(fc_points, field_discharge_mean, 'DOUBLE')
    arcpy.AddField_management(fc_points, field_discharge_min, 'DOUBLE')
    arcpy.AddField_management(fc_points, field_head, 'FLOAT')
    arcpy.AddField_management(fc_points, field_power, 'DOUBLE')
    arcpy.AddField_management(fc_points, field_power_max, 'DOUBLE')
    arcpy.AddField_management(fc_points, field_power_mean, 'DOUBLE')
    arcpy.AddField_management(fc_points, field_power_min, 'DOUBLE')
    with arcpy.da.UpdateCursor(fc_points, [field_discharge_accumulated, field_discharge_max, field_discharge_mean,
                                           field_discharge_min, field_head, field_power, field_power_max,
                                           field_power_mean, field_power_min]) as cursor:
        for index, row in enumerate(cursor, start=0):
            row[0] = points_df.loc[index, field_discharge_accumulated]
            row[1] = points_df.loc[index, field_discharge_max]
            row[2] = points_df.loc[index, field_discharge_mean]
            row[3] = points_df.loc[index, field_discharge_min]
            row[4] = points_df.loc[index, field_head]
            row[5] = points_df.loc[index, field_power]
            row[6] = points_df.loc[index, field_power_max]
            row[7] = points_df.loc[index, field_power_mean]
            row[8] = points_df.loc[index, field_power_min]
            cursor.updateRow(row)
            print('added {}'.format(index))


def runner():
    arcpy.env.workspace = r'C:/Users/Chris/Documents/GIS/HydroModeller.gdb'
    fc_rivers = 'rivers_ug_small_projected'
    fc_flow_acc = 'flow_acc_af_projected'
    fc_gscd = 'gscd_qmean_projected'
    fc_land_type = 'land_type_projected'
    # TODO maybe need to aggregate elevation, as doesn't capture always
    fc_elevation = 'elevation_af_projected'
    fc_points = 'hydro_points_ug'
    fc_et_ref = 'et_ref_projected'

    arcpy.env.overwriteOutput = True
    arcpy.env.addOutputsToMap = False
    interval = 1000

    network, nodes = create_network(fc_rivers)
    network = calc_stream_order(network, nodes)
    nodes_df = prepare_nodes(fc_rivers, nodes, fc_flow_acc, fc_gscd, fc_land_type, fc_et_ref)
    nodes_df = calc_flow_acc(nodes_df, network)
    nodes_df, runoff_df = rainfall_runoff(nodes_df)
    network_df, nodes_df, discharge_df = calc_discharge(nodes_df, runoff_df, network, fc_flow_acc, fc_rivers)
    calc_hydro_potential(fc_rivers, fc_elevation, interval, fc_points, network_df)


if __name__ == "__main__":
    startTime = datetime.now()
    runner()
    print('time {}'.format(datetime.now() - startTime))
