from __future__ import absolute_import, division, print_function, unicode_literals
import arcpy
import numpy as np
from datetime import datetime


def shreve(arc_index, direction_node_id, network, nodes):
    """
    Caclulate Shreve stream order instead of Strahler
    This ensures that a downstream river is always higher
    """
    import heapq
    up_stream_orders = []
    if len(nodes[direction_node_id]) == 4:
        network[arc_index][3] = 1
    else:
        for index, arc in enumerate(nodes[direction_node_id]):
            if index >= 3:
                if network[arc][0] != arc_index:
                    if network[arc][1] != direction_node_id:
                        up_stream_orders.append(shreve(arc, network[arc][1], network, nodes))
                    else:
                        up_stream_orders.append(shreve(arc, network[arc][2], network, nodes))

        max_orders = heapq.nlargest(2, up_stream_orders)
        order = 0 + max_orders[0] + max_orders[1]

        network[arc_index][3] = order

    print('so' + str(arc_index))
    return network[arc_index][3]


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
        network[arc_index][3] = 1
    else:
        for arc in nodes[direction_node_id]:
            if network[arc][0] != arc_index:
                if network[arc][1] != direction_node_id:
                    up_stream_orders.append(strahler(arc, network[arc][1], network, nodes))
                else:
                    up_stream_orders.append(strahler(arc, network[arc][2], network, nodes))
        max_order = 0
        max_order_count = 0
        for order in up_stream_orders:
            if order > max_order:
                max_order = order
                max_order_count = 1
            elif order == max_order:
                max_order_count += 1
        if max_order_count > 1:
            network[arc_index][3] = max_order + 1
        else:
            network[arc_index][3] = max_order
    print('so' + str(arc_index))
    return network[arc_index][3]


def create_network(rivers):
    # import the start and end end points into the network array of arcs
    # and simultaneously import them into a separate node array (so two nodes per arc)
    # there will be duplicate nodes, but they duplicates are ignored later
    counter = 0
    length = int(arcpy.GetCount_management(rivers).getOutput(0))
    network = np.empty([length, 6], dtype=int)
    nodes = []
    with arcpy.da.SearchCursor(rivers, ['SHAPE@']) as cursor:
        for row in cursor:
            # the geometry variables
            fx = int(row[0].firstPoint.X)
            fy = int(row[0].firstPoint.Y)
            lx = int(row[0].lastPoint.X)
            ly = int(row[0].lastPoint.Y)

            # store the index as a column for easier access
            # the last column is for stream order
            network[counter] = [counter, fx, fy, lx, ly, 0]

            # create nodes for the start point and end point separately
            nodes.append([2*counter, fx, fy])
            nodes.append([2*counter+1, lx, ly])

            counter += 1
            print('imported ' + str(counter))

    # the most time consuming part of the script
    # runs through every arc and every node and where there is a match, add
    #  a reference to the index into the network array
    # once both nodes have been found for an arc, break that for loop and
    #  start with the next arc
    for arc in network:
        breaker = 0
        for node in nodes:
            # if the fx and fy values match
            if arc[1] == node[1] and arc[2] == node[2]:
                # overwrite fx with an index
                arc[1] = node[0]
                breaker += 1
            # if the lx and ly values match
            elif arc[3] == node[1] and arc[4] == node[2]:
                # overwrite lx with an index
                arc[3] = node[0]
                breaker += 1
            if breaker > 1:
                # reset and skip to the next arc
                break
        print('indexed ' + str(arc[0]))

    # delete all the X and Y values from network array to make next steps simpler
    network = np.delete(network, [2, 4], 1)
    print('deleted stuff')

    # for every node, add references to every arc that connects to it
    for arc in network:
        # tell the arc's starting node that it exists
        nodes[arc[1]].append(arc[0])
        # and tell the arc's ending node that it exists
        nodes[arc[2]].append(arc[0])
        print('nodes ' + str(arc[0]))

    return network, nodes, length


def get_discharge(rivers, network, nodes, length, flowacc_fc, gscd_fc):
    # Create a points feature class, and add points using the coordinates of every node
    desc = arcpy.Describe(rivers)
    nodes_fc = rivers + '_temp_nodes'
    arcpy.CreateFeatureclass_management(arcpy.env.workspace, nodes_fc, geometry_type='POINT',
                                        spatial_reference=desc.spatialReference)
    with arcpy.da.InsertCursor(nodes_fc, ['SHAPE@XY']) as cursor:
        for node in nodes:
            xy = (node[1], node[2])
            cursor.insertRow([xy])

    # Extract the Flow accumulation and gscd values into these points, from the rasters
    flowacc_field = 'FlowAccummulation'
    gscd_field = 'RunOff'

    arcpy.CheckOutExtension("Spatial")

    # Created an upsampled version of the flowacc feature class, to make sure nodes hit the value
    flowacc_agg = arcpy.sa.Aggregate(flowacc_fc, 10, "MAXIMUM", "EXPAND", "DATA")
    flowacc_agg_fc = 'flowacc_agg_temp'
    flowacc_agg.save(flowacc_agg_fc)

    # might have to run these directly in ArcMap (license issue)
    # sa.ExtractMultiValuesToPoints(nodes_fc, [[flowacc_fc, flowacc_field]])
    arcpy.sa.ExtractMultiValuesToPoints(nodes_fc, [[flowacc_agg_fc, flowacc_field]])
    arcpy.sa.ExtractMultiValuesToPoints(nodes_fc, [[gscd_fc, gscd_field]])

    # Create an array called node_dis
    # It has the following structure:
    # - 0       - 1         - 2         - 3                 - 4         -
    # - index   - flowacc   - runoff    - actual_flowacc    - discharge -
    # load flowacc and runoff from the shapefile we just created
    node_dis = np.zeros([len(nodes), 5], dtype=float)
    counter = 0
    with arcpy.da.SearchCursor(nodes_fc, [flowacc_field, gscd_field]) as cursor:
        for row in cursor:
            node_dis[counter] = [counter, row[0], row[1], 0, 0]
            print('added nodedis ' + str(counter))
            counter += 1

    # Clean up the coordinates from nodes now that we're done with that
    nodes = [row[3:] for row in nodes]

    # For each node, subtract the flow accumulation for all upstream nodes,
    # so that we calculate it's 'self discharge' only for it's directly contributing area
    for index, node in enumerate(nodes):
        actual_flowacc = node_dis[index][1]  # flowAcc
        for arc in node:
            if network[arc][2] == index:  # up the arc flows *into* the node
                subtract = node_dis[network[arc][1]][1]
                actual_flowacc -= subtract  # subtract upstream flowAcc
        node_dis[index][3] = actual_flowacc

    cell = arcpy.GetRasterProperties_management(flowacc_fc, 'CELLSIZEX').getOutput(0)
    area = float(cell) * float(cell)
    # and calculate the proper value in m^3/s - 8760 (h/year) * 3600 (s/h) * 1000 (mm/m) = 3153600000

    # calculate self discharge for each node
    # some flowacc values are negative, because hydrosheds doesn't line up properly
    # therefore, any points that are negative set to use their original flowaccu
    for node in node_dis:
        node[4] = node[2]*node[3]*area/3153600000
        if node[4] < 0:
            node[4] = 0

    # start from Shreve order 1, and contribute all discharges from upstread nodes to downstream nodes
    for so in range(1, max(network.T[3])+1):
        for arc in network:
            if arc[3] == so:
                node_dis[arc[2]][4] += node_dis[arc[1]][4]*0.8

    # create a new array indexed to the network array
    # and add the discharge from the upstream node of each arc to that arc
    discharge = np.empty(length, dtype=float)
    for index, arc in enumerate(network):
        discharge[index] = node_dis[arc[1]][4]

    # add the results back into the ArcGIS feature
    discharge_field = 'Discharge'
    arcpy.AddField_management(rivers, discharge_field, 'FLOAT')
    with arcpy.da.UpdateCursor(rivers, [discharge_field]) as cursor:
        for index, row in enumerate(cursor, start=0):
            row[0] = discharge[index]
            cursor.updateRow(row)
            print('added ' + str(index))


def calc_hydro_potential(rivers, elevation, discharge, interval, points):
    desc = arcpy.Describe(rivers)
    arcpy.CreateFeatureclass_management(arcpy.env.workspace, points, geometry_type='POINT',
                                        spatial_reference=desc.spatialReference)
    fid_name = 'river_ID'
    arcpy.AddField_management(points, fid_name, 'LONG')

    with arcpy.da.SearchCursor(rivers, ['SHAPE@', desc.OIDFieldName]) as search_cursor:
        with arcpy.da.InsertCursor(points, ['SHAPE@', fid_name]) as insert_cursor:
            for row in search_cursor:
                line = row[0]

                if line:
                    cur_length = interval
                    max_position = line.length
                    insert_cursor.insertRow([line.firstPoint, row[1]])
                    while cur_length < max_position:
                        insert_cursor.insertRow([line.positionAlongLine(cur_length, False), row[1]])
                        cur_length += interval

    arcpy.sa.ExtractMultiValuesToPoints(points, [[elevation, 'elev']])

    # calculate the elevation difference between successive points
    arcpy.AddField_management(points, 'head', 'FLOAT')
    with arcpy.da.UpdateCursor(points, [fid_name, 'elev', 'head']) as cursor:
        prev_fid = 0
        prev_elev = 0

        for row in cursor:
            if row[1]:  # check to make sure an elevation entry exists
                current_fid = row[0]
                current_elev = row[1]

                if (current_fid == prev_fid) and ((prev_elev - current_elev) > 0):
                    row[2] = prev_elev - current_elev
                else:
                    row[2] = 0

                cursor.updateRow(row)
                prev_fid = current_fid
                prev_elev = current_elev
            else:
                row[2] = 0
                prev_fid = row[0]
                prev_elev = 0

    # delete all points where head is 0 or negative
    with arcpy.da.UpdateCursor(points, 'head') as cursor:
        for row in cursor:
            if row[0] <= 0:
                cursor.deleteRow()

    arcpy.sa.ExtractMultiValuesToPoints(points, [[discharge, 'discharge_second']])

    # before calculating the power, delete any points that will result in a value absurdly large or small
    # (LONG only goes up to 2 billion)
    # currently set at less than 10 litres/sec, or more than 2000 m^3/sec
    discharge_lower_cutoff = 0.01
    discharge_upper_cutoff = 2000
    with arcpy.da.UpdateCursor(points, 'discharge_second') as cursor:
        for row in cursor:
            if row[0] < discharge_lower_cutoff or row[0] > discharge_upper_cutoff:
                cursor.deleteRow()

    # calculate the power in watts based on Alex's formula
    # P = rho * g * nt * ng * conv * Q * deltaH
    rho = '1000'  # density
    g = '9.81'  # gravity
    nt = '0.88'  # turbine efficiency
    ng = '0.96'  # generator efficiency
    conv = '0.6'  # conversion factor for environmental flow deduction
    arcpy.AddField_management(points, 'power', 'FLOAT')
    arcpy.CalculateField_management(points, 'power',
                                    rho + '*' + g + '*' + nt + '*' + ng + '*' + conv + '* !discharge_second! * !head!',
                                    'PYTHON_9.3')


def runner():
    arcpy.env.workspace = r'C:/Users/Chris/Work/GIS/Data/Default.gdb'
    rivers = 'riv_beni'
    flowacc_fc = r'C:/Users/Chris/Work/GIS/Data/Hydrology/elevation/flow_accu'
    gscd_fc = r'C:/Users/Chris/Work/GIS/Data/Hydrology/streamflow/gscd'

    arcpy.env.overwriteOutput = True
    arcpy.env.addOutputsToMap = False
    object_id = 5
    sink = object_id - 1

    network, nodes, length = create_network(rivers)

    # run the stream ordering function
    # sink is the index of the sink arc
    network[sink][3] = shreve(sink, network[sink][1], network, nodes)

    get_discharge(rivers, network, nodes, length, flowacc_fc, gscd_fc)

if __name__ == "__main__":
    startTime = datetime.now()
    runner()
    print('time ' + str(datetime.now() - startTime))
