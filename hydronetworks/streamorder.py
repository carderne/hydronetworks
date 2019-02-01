# streamorder.py
#! python3

"""Functions for calculating stream order."""

import heapq


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


def shreve(arc_index, direction_node_id, network, nodes):
    """
    Caclulate Shreve stream order instead of Strahler
    This ensures that a downstream river is always higher
    """

    up_stream_orders = []
    if len(nodes[direction_node_id]) == 5:
        network[arc_index][7] = 1
    else:
        for index, arc in enumerate(nodes[direction_node_id]):
            if index >= 4:
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

    # print('so {}'.format(arc_index))
    return network[arc_index][7]

