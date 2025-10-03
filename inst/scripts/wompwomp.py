import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata
import igraph
import re
from pyalluvial.alluvial import sigmoid, calc_sigmoid_line
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import itertools
import warnings
from python_tsp.exact import solve_tsp_dynamic_programming
from sklearn.metrics import adjusted_rand_score
from typing import Tuple, List
from splitspy.nnet import nnet_cycle, nnet_splits
import random

# function to plot the alluvium (adapted from pyalluvial)
def _plot_alluvium(df, xaxis_names, y_name='value', alluvium=None, 
                   color_alluvium = False, color_alluvium_boundary = False,
                   order_dict=None, all_color_dict=None, ax=None, x_init=0,
                  objective_calc=False, invert_xy = False, alluvial_alpha = 0.5,
                  box_line_width = 1,box_width = 0.4, alluvial_edge_width = 0.1):
    # if no specified reference, make column 0 the alluvial reference
    if alluvium is None:
        alluvium = xaxis_names[0]
    else:
        color_alluvium = True
        
    if ax is None and not objective_calc:
        warnings.warn("No ax object but not performing objective calculations")
        objective_calc == True
        
    if order_dict is None:
        pass
    else:
        for key, orders in order_dict.items():
            if key in df.columns:
                # reverse order_dict contents to plot as is intuitive
                df[f'{key}_order'] = df[key].map(dict(zip(orders[::-1], list(range(len(orders)))))).astype(int)
                #df[f'{key}_order'] = df[key].map(dict(zip(orders, list(range(len(orders)))))).astype(int)
    xaxis_after_ordering = []
    for xaxis in xaxis_names:
        if xaxis + '_order' in df.columns:
            xaxis_after_ordering += [xaxis+ '_order']
        else:
            xaxis_after_ordering += [xaxis]
    df = df.sort_values(xaxis_after_ordering).reset_index(drop=True).reset_index(names='alluvia')

    for xaxis in xaxis_names:
        if xaxis + '_order' in df.columns:
            corrected_name = xaxis + '_order'
            sorting_list = [corrected_name] + [x for x in xaxis_after_ordering if x != corrected_name]
            df[f'y_{xaxis}'] = df.sort_values(sorting_list)[y_name].cumsum().shift(1).fillna(0)
        else:
            sorting_list = [xaxis] + [x for x in xaxis_after_ordering if x != xaxis]
            df[f'y_{xaxis}'] = df.sort_values(sorting_list)[y_name].cumsum().shift(1).fillna(0)
    if objective_calc:
        return df
    
    # plot alluvium
    for i in range(len(xaxis_names)):
        if i + 1 >= len(xaxis_names):
            break
        else:
            for y_left, y_right, height, color_key in zip(
                df[f'y_{xaxis_names[i]}'].values,
                df[f'y_{xaxis_names[i + 1]}'].values,
                df[y_name].values,
                df[alluvium].values):
                
                xs, ys_under = calc_sigmoid_line(1 - box_width, y_left, y_right)
                xs += i + box_width/2
                ys_upper = ys_under + height
                if not objective_calc:
                    label_key = color_key
                    if label_key in ax.get_legend_handles_labels()[1]:
                        label_key = None
                    if invert_xy:
                        if all_color_dict is None or not color_alluvium:
                            ax.fill_betweenx(xs + x_init, ys_under, ys_upper, alpha=alluvial_alpha, color='grey', 
                                             edgecolor='black', label = label_key,
                                                linewidth = alluvial_edge_width)
                        else:
                            ax.fill_betweenx(xs + x_init, ys_under, ys_upper, alpha=alluvial_alpha, 
                                                 color=all_color_dict[color_key], 
                                                 edgecolor=all_color_dict[color_key] if color_alluvium_boundary else 'black', 
                                                 label = label_key,
                                                linewidth = alluvial_edge_width)
                    else:
                        if all_color_dict is None or not color_alluvium:
                            ax.fill_between(xs + x_init, ys_under, ys_upper, alpha=alluvial_alpha, color='grey', 
                                            edgecolor='black', label = label_key,
                                                linewidth = alluvial_edge_width)
                        else:
                            ax.fill_between(xs + x_init, ys_under, ys_upper, alpha=alluvial_alpha, 
                                                 color=all_color_dict[color_key], 
                                                 edgecolor=all_color_dict[color_key] if color_alluvium_boundary else 'black', 
                                                 label = label_key,
                                                linewidth = alluvial_edge_width)



def plot_alluvial_internal(df, xaxis_names, y_name, alluvium_column = None, order_dict=None, color_dict = None,
         color_boxes = True, color_alluvium = False, alluvial_alpha = 0.5, color_alluvium_boundary = False,
         ignore_continuity=False, cmap_name='tab20', 
        include_labels_in_boxes = True, min_text = 4, default_text_size = 14, autofit_text = True,
        box_line_width = 1,box_width = 0.4, alluvial_edge_width = 0.1,
         figsize=(6.4, 4.8),
        y_axis_label = False, invert_xy = False, include_stratum_legend = False, include_alluvium_legend = False, verbose = False
        ):
    '''
    Plot alluvial plot
    
    Parameters
    --------
    df: pandas DataFrame
        df must be in wide format, as shown in the example below.
        Example)
        |    |   survived | class   | sex    |   freq |
        |---:|-----------:|:--------|:-------|-------:|
        |  0 |          0 | First   | female |      3 |
        |  1 |          0 | First   | male   |     77 |
        |  2 |          0 | Second  | female |      6 |
        |  3 |          0 | Second  | male   |     91 |
        |  4 |          0 | Third   | female |     72 |
    xaxis_names: list of str
        Specify the column names to line up on the x-axis. It is drawn in the order of this list.
        Example) ['class', 'sex']
    y_name: str
        Specify a column name that represents the height of the y-axis.
        Example) 'freq'
    alluvium_column: str or None
        Specify the column name of the alluvial color. 
        If ignore_continuity is true and you want to reset the colors in each Stratum, set it to None.
        Example) 'survived'
    order_dict: dict
        If you want to adjust the order in each Stratum, specify the order like the following example.
        Example) {'class': ['Third', 'Second', 'First'], 'sex': ['male', 'female']}
    ignore_continuity: bool
        Specify True if you want to ignore the continuity of each axis, otherwise False.
        Example) True
    color_dict: dict
        manual selection of colors 
    cmap_name: str
        Specify any matplotlib's colormap name you want to use in the following link.
        It is recommended to choose from the Qualitative colormaps.
        ref) https://matplotlib.org/examples/color/colormaps_reference.html
        Example) 'tab10'
    figsize: tuple of float
        Specify the figsize.
        Example) (10, 10)
        
    Return
    --------
    fig: matplotlib figure object
    '''
    df = df.copy()
    if verbose:
        print('Assigning colors')
    if color_boxes or alluvium_column is not None:
        if color_dict is None:
            color_dict = {}
        if alluvium_column in xaxis_names or alluvium_column is None:
            keys = pd.melt(df[xaxis_names])['value'].unique().tolist()
        else:
            keys = pd.melt(df[xaxis_names+ [alluvium_column]])['value'].unique().tolist()
        if len(color_dict.keys()) < len(keys):
            cmap = plt.get_cmap(name=cmap_name, lut=len(keys))
            available_colors = cmap(np.arange(0,cmap.N)) 
            for index, key in enumerate(keys):
                if key not in color_dict.keys():
                    color_dict[key] = available_colors[index]
    if verbose:
        print('Plotting Stratum')
        
    # 各stratumの高さを計算する
    stratum_dict = {}
    for xaxis in xaxis_names:
        stratum_dict[xaxis] = df.groupby(xaxis)[y_name].sum()
    
    # stratumの順番を設定する
    if order_dict is None:
        pass
    else:
        for key, orders in order_dict.items():
            stratum = stratum_dict[key]
            # reverse order_dict contents to plot as is intuitive
            stratum_dict[key] = stratum[orders[::-1]]
            #stratum_dict[key] = stratum[orders]

    fig, ax = plt.subplots(figsize=figsize)

    stratum_labels = {}
    # plot stratum (stacked bar)
    for i, stratum in enumerate(stratum_dict.values()):
        xtick_label = stratum.index.name
        names = stratum.index.values
        values = stratum.values
        stratum_labels[i] = {}

        for j, (name, value) in enumerate(zip(names, values)):
            bottom = values[:j - len(stratum)].sum()  
            if invert_xy:     
                rectangle = ax.barh(
                        y=[i],
                        width=value,
                        left=bottom,
                        color=color_dict[name] if color_boxes else 'white',
                        label = name if include_stratum_legend else None,
                        edgecolor='black',
                        fill=True,
                        linewidth=box_line_width,
                        height=box_width
                    )
                text_y = rectangle[0].get_width() / 2 + bottom
                
            else:            
                rectangle = ax.bar(
                        x=[i],
                        height=value,
                        bottom=bottom,
                        color=color_dict[name] if color_boxes else 'white',
                        label = name if include_stratum_legend else None,
                        edgecolor='black',
                        fill=True,
                        linewidth=box_line_width,
                        width=box_width
                    )
                text_y = rectangle[0].get_height() / 2 + bottom
            stratum_labels[i][name] = [rectangle[0], text_y]
    
    
    if verbose:
        print('Plotting Alluvia')
    if ignore_continuity:
        for i in range(len(xaxis_names)):
            if i + 1 >= len(xaxis_names):
                break
            else:
                if alluvium_column is None:
                    alluvium_column = xaxis_names[i]
                    color_val = df[alluvium_column].unique()
                    color_dict = dict(zip(color_val, list(range(len(color_val)))))
                
                agg_cols = list(set([f'{alluvium}'] + xaxis_names[i: i+2]))
                df_agg = df.groupby(agg_cols, as_index=False)[y_name].sum()

                
                _plot_alluvium(df_agg, xaxis_names[i: i+2], y_name, alluvium = alluvium_column, 
                               color_alluvium = color_alluvium, color_alluvium_boundary = color_alluvium_boundary,
                               order_dict=order_dict, all_color_dict = color_dict, ax=ax, x_init=i,
                              invert_xy=invert_xy, alluvial_alpha = alluvial_alpha, alluvial_edge_width = alluvial_edge_width)
                alluvium = None
    else:
        _plot_alluvium(df, xaxis_names, y_name, alluvium = alluvium_column, 
                       color_alluvium = color_alluvium, color_alluvium_boundary = color_alluvium_boundary,
                       order_dict=order_dict, all_color_dict = color_dict, ax=ax, x_init=0,
                      invert_xy=invert_xy, alluvial_alpha = alluvial_alpha,
                      box_line_width = box_line_width,box_width = box_width, alluvial_edge_width = alluvial_edge_width)
    
    #ax.set_axis_off()
    if invert_xy:
        if not y_axis_label:
            ax.get_xaxis().set_visible(False)    
            ax.spines['bottom'].set_visible(False)  
        # set xticklabels
        ax.set_yticks(list(range(len(stratum_dict))))
        ax.set_yticklabels([stratum.index.name for stratum in stratum_dict.values()])
    else:
        if not y_axis_label:
            ax.get_yaxis().set_visible(False)    
            ax.spines['left'].set_visible(False)  
        # set xticklabels
        ax.set_xticks(list(range(len(stratum_dict))))
        ax.set_xticklabels([stratum.index.name for stratum in stratum_dict.values()])


    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    if verbose:
        print('Adding Stratum labels')
    # add box labels
    if include_labels_in_boxes:
        for xaxis in stratum_labels.keys():
            for name in stratum_labels[xaxis]:
                rect, y = stratum_labels[xaxis][name]
                bar_label = ax.text(x=xaxis, y=y, s=name, horizontalalignment='center', verticalalignment='center',
                                    wrap=True, size = default_text_size)
                bar_label._get_wrap_line_width = lambda : rect.get_width()
                if autofit_text:
                    auto_fit_fontsize(text=bar_label, rectangle=rect, fig=fig, ax=ax, text_min = min_text)
    if include_stratum_legend or include_alluvium_legend:
        lgd = plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        place_legend(fig, ax, lgd)
    return fig, ax

def place_legend(fig, ax, lgd):
    plt.gcf().draw_without_rendering()
    invFigure = plt.gcf().transFigure.inverted()
    lgd_pos = lgd.get_window_extent()
    lgd_coord = invFigure.transform(lgd_pos)
    lgd_xmax = lgd_coord[1, 0]
    ax_pos = plt.gca().get_window_extent()
    ax_coord = invFigure.transform(ax_pos)
    ax_xmax = ax_coord[1, 0]
    shift = 1 - (lgd_xmax - ax_xmax)
    plt.gcf().tight_layout(rect=(0, 0, shift, 1))

def auto_fit_fontsize(text,rectangle, fig, ax, text_min=1):
    '''Auto-decrease the fontsize of a text object.

    Args:
        text (matplotlib.text.Text)
        width (float): allowed width in data coordinates
        height (float): allowed height in data coordinates
    '''# get text bounding box in figure coordinates
    renderer = fig.canvas.get_renderer()
    fig.draw_without_rendering()
    
    bbox_text = text.get_window_extent(renderer=renderer)
    rect_bbox_text = rectangle.get_window_extent(renderer=renderer)
    # transform bounding box to data coordinates
    bbox_text = Bbox(ax.transAxes.inverted().transform(bbox_text))
    rect_bbox_text = Bbox(ax.transAxes.inverted().transform(rect_bbox_text))
    
    # evaluate fit and recursively decrease fontsize until text fits
    fits_width = bbox_text.width < rect_bbox_text.width 
    fits_height = bbox_text.height < rect_bbox_text.height 
    
    if not all((fits_width, fits_height)):
        text.set_fontsize(text.get_fontsize()-1)
        if text.get_fontsize() > text_min:
            auto_fit_fontsize(text, rectangle, fig, ax)


def __setup_matrix_custom(labels: [str], matrix: [float]) -> np.array:
    n = len(labels)
    max_number_of_nodes = max(3, 3 * n - 5)

    # mat = np.empty(((max_number_of_nodes + 1), (max_number_of_nodes + 1)))  # JMR, commented out (to remove randomness)
    mat = np.zeros((max_number_of_nodes + 1, max_number_of_nodes + 1), dtype=np.float64)  # JMR, commented out (to remove randomness)

    for i in range(0, max_number_of_nodes):
        for j in range(0, max_number_of_nodes):
            if i < n and j < n:
                mat[i+1][j+1] = matrix[i][j]

    return mat

def compute_custom(labels: [str], matrix) -> [int]:
    n = len(labels)

    if n <= 3:
        return list(range(0, n + 1))

    nodes_head = nnet_cycle.__setup_nodes(n)

    mat = __setup_matrix_custom(labels, matrix.copy())  # matrix is 0-based, mat is 1-based

    joins = nnet_cycle.__join_nodes(n, mat, nodes_head)

    cycle = nnet_cycle.__expand_nodes(joins, nodes_head)

    cycle = nnet_cycle.__normalize_cycle(cycle)

    return cycle

def neighbor_net(labels: List[str], mat: List[List[float]], cutoff=0.0001, constrained=True, compute_splits=False):
    mat = mat.astype(np.float64)
    cycle = compute_custom(labels, mat)
    if compute_splits:
        splits = nnet_splits.compute(len(labels), mat, cycle, cutoff, constrained)
    else:
        splits = []
    # correct for 1-indexing in nnet
    cycle = [x-1 for x in cycle[1:]]
    return cycle, splits


def generate_distance_matrix(df, graphing_columns,
                    col_weights = 'value', matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, 
                     weight_scalar = 5e5, nn_normalize = True):

    matrix_df = df[graphing_columns].astype(str).apply(lambda x : x.name+'~~'+x)
    matrix_df[col_weights] = df[col_weights]
    all_nodes = pd.melt(matrix_df[graphing_columns])['value'].unique().tolist()
    
    full_dist_matrix = np.full((len(all_nodes),len(all_nodes)), matrix_initialization_value)
    
    if matrix_initialization_value != same_side_matrix_initialization_value:
        for col in graphing_columns:
            col_indices = list(np.where([col in x for x in all_nodes])[0])
            for x in col_indices:
                full_dist_matrix[x, col_indices] = same_side_matrix_initialization_value
                full_dist_matrix[col_indices, x] = same_side_matrix_initialization_value
    pairwise_groupings = list(itertools.combinations(graphing_columns, 2))
    
    results = []
    
    for cols in pairwise_groupings:
        grouped = df.groupby(list(cols), as_index=False).agg(total_value=(col_weights, 'sum'))
        grouped['grouping'] = '+'.join(list(cols))
        results.append(grouped)
    final_result = pd.concat(results).reset_index(drop=True)
    final_result = final_result[final_result['total_value'] > 0]
    final_result[['comp1', 'comp2']] = final_result.grouping.str.split('+', expand=True)
    final_result['comp1_results'] = final_result.values[final_result.index.get_indexer(final_result['comp1'].index), final_result.columns.get_indexer(final_result['comp1'])]
    final_result['comp2_results'] = final_result.values[final_result.index.get_indexer(final_result['comp2'].index), final_result.columns.get_indexer(final_result['comp2'])]
    
    final_result['comp1'] = final_result['comp1'].astype('str') + '~~' + final_result['comp1_results'].astype('str')
    final_result['comp2'] = final_result['comp2'].astype('str') + '~~' + final_result['comp2_results'].astype('str')
    
    final_result['loc1'] = final_result.comp1.apply(lambda x: all_nodes.index(x))
    final_result['loc2'] = final_result.comp2.apply(lambda x: all_nodes.index(x))
    
    for index, row in final_result.iterrows():
            if nn_normalize:
                full_dist_matrix[final_result['loc1'], final_result['loc2']] = weight_scalar * -np.log(final_result['total_value'])
                full_dist_matrix[final_result['loc2'], final_result['loc1']] = weight_scalar * -np.log(final_result['total_value'])
            else:
                full_dist_matrix[final_result['loc1'], final_result['loc2']] = weight_scalar * final_result['total_value']
                full_dist_matrix[final_result['loc2'], final_result['loc1']] = weight_scalar * final_result['total_value']
    
    if nn_normalize:
         # make sure all numbers are positive for neighbornet
        min_val_abs = np.abs(np.min(full_dist_matrix))
        full_dist_matrix = full_dist_matrix + (min_val_abs + 1)
        full_dist_matrix[np.isnan(full_dist_matrix)] = 1e6
    
    return full_dist_matrix, all_nodes

def sort_dist_matrix(full_dist_matrix, all_nodes, sorting_algorithm = "neighbornet"):
    if sorting_algorithm == 'neighbornet':
        result = neighbor_net(all_nodes, full_dist_matrix) 
        cycle = pd.Series(all_nodes)[result[0]].tolist()
    elif sorting_algorithm == 'tsp':
        permutation, distance = solve_tsp_dynamic_programming(full_dist_matrix)
        cycle = pd.Series(all_nodes)[permutation].tolist()
    else:
        cycle = all_nodes
    return cycle

def rotate_left(vec, k=1):
    n = len(vec)
    k = k%n
    if k !=0:
        vec = vec[k:n] + vec[0:k]
    return vec

def get_order_dict(cycle, graphing_columns=None):
    graph_df = pd.Series(cycle).str.split('~~', expand = True)
    graphs = {}
    if graphing_columns is None:
        for key in graph_df[0].unique():
            graphs[key] = graph_df[graph_df[0] == key][1].values.tolist()
    else:
        for key in graphing_columns:
            graphs[key] = graph_df[graph_df[0] == key][1].values.tolist()
    return graphs

# calculate_objective.py
class BIT:
    def __init__(self, size):
        self.tree_count = [0] * (size + 1)
        self.tree_weight = [0.0] * (size + 1)

    def update(self, index, weight):
        index += 1
        while index < len(self.tree_count):
            self.tree_count[index] += 1
            self.tree_weight[index] += weight
            index += index & -index

    def query(self, index):
        count = 0
        weight_sum = 0.0
        index += 1
        while index > 0:
            count += self.tree_count[index]
            weight_sum += self.tree_weight[index]
            index -= index & -index
        return count, weight_sum

    def query_range(self, low, high):
        c1, w1 = self.query(high)
        c2, w2 = self.query(low - 1)
        return c1 - c2, w1 - w2

def calculate_objective_fenwick(df, y1, y2, weight = 'count'):
    # Step 1: Sort by y1
    df_sorted = df.sort_values(y1).reset_index(drop=True)

    # Step 2: Rank-compress y2 (higher y2 → higher rank)
    df_sorted['y2_rank'] = rankdata(df_sorted[y2], method='dense').astype(int)
    max_rank = df_sorted['y2_rank'].max()

    # Step 3: Use BIT
    bit = BIT(size=max_rank + 2)
    total_cross_weight = 0.0

    for y2_rank, weight in zip(df_sorted['y2_rank'].values, df_sorted[weight].values):
        # Count previous y2s > current (strictly greater)
        _, sum_weights_above = bit.query_range(y2_rank + 1, max_rank)
        total_cross_weight += weight * sum_weights_above

        # Add current y2_rank to BIT
        bit.update(y2_rank, weight)
    
    return total_cross_weight

def determine_crossing_edges(df, graphing_columns = None, order_dict = None,
                             col_weights = "value"):
    objective_val = 0
    lode_df = _plot_alluvium(df, graphing_columns, col_weights, order_dict=order_dict,
                      objective_calc=True)
    lode_df = lode_df.sort_index()
    lode_df = lode_df.rename_axis('alluvium').reset_index()
    for h in range(len(graphing_columns)-1):
        x1 = graphing_columns[h]
        x2 = graphing_columns[h+1]
    
        y1 = 'y_' + graphing_columns[h]
        y2 = 'y_' + graphing_columns[h+1]
        
        objective_val += calculate_objective_fenwick(lode_df, y1, y2, weight=col_weights)
    return objective_val

def determine_column_order(df_gather, cycle, graphing_columns, column_weights = 'value',
                           matrix_initialization_value_column_order = 1e6, 
                        weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", 
                        column_sorting_algorithm = "tsp", 
                                 verbose = False):
    full_dist_matrix = np.full((len(graphing_columns),len(graphing_columns)), matrix_initialization_value_column_order)
    
    pairwise_groupings = list(itertools.combinations(graphing_columns, 2))
    for i in range(len(pairwise_groupings)):
        col1 = pairwise_groupings[i][0]
        col2 = pairwise_groupings[i][1]

        col1_index = np.where(pd.Series(graphing_columns) == col1)[0][0]
        col2_index = np.where(pd.Series(graphing_columns) == col2)[0][0]

        if column_sorting_metric == "edge_crossing":
            neighbornet_objective = determine_crossing_edges(df_gather, graphing_columns = [col1, col2], 
                                     order_dict = get_order_dict(cycle, [col1, col2]), col_weights = column_weights)
            neighbornet_objective = weight_scalar_column_order * np.log1p(neighbornet_objective)
        elif column_sorting_metric == "ARI":
            expanded_df = df_gather.loc[df_gather.index.repeat(df_gather[column_weights])]
            neighbornet_objective = adjusted_rand_score(expanded_df[col1], expanded_df[col2])
            neighbornet_objective = 1 - neighbornet_objective
            neighbornet_objective = weight_scalar_column_order * 50 * neighbornet_objective
        else:
            ValueError(f"{column_sorting_metric} is not a valid option")
            
        full_dist_matrix[col1_index, col2_index] = neighbornet_objective
        full_dist_matrix[col2_index, col1_index] = neighbornet_objective

    nodes = graphing_columns
    result = sort_dist_matrix(full_dist_matrix, nodes, sorting_algorithm = column_sorting_algorithm)
    n = len(result)
    adj_distances = [full_dist_matrix[x, x%(n-1)+1] for x in range(n)]
    result = rotate_left(result, np.argmax(adj_distances)+1)
    
    return result


def determine_optimal_cycle_start(df_gather, cycle, graphing_columns, column_weights = 'value', optimize_column_order = True,
                                 optimize_column_order_per_cycle = False, cycle_start_positions = None,
                                  matrix_initialization_value_column_order = 1e6, 
                        weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", 
                        column_sorting_algorithm = "neighbornet", 
                                 verbose = False):
    neighbornet_objective_minimum = np.inf
    cycle_best = cycle
    
    for i in range(len(cycle)):
        if (cycle_start_positions is not None and (i + 1) not in cycle_start_positions):
            continue
            
        # Rotate stratum order and determine neighbornet objective
        cycle_shifted = rotate_left(cycle, i)
        neighbornet_objective = determine_crossing_edges(df_gather, graphing_columns, get_order_dict(cycle_shifted), col_weights=column_weights)
        
        if verbose:
            print(f"neighbornet_objective for iteration {i} = {neighbornet_objective}")
            
        # if better, update best cycle and minimum
        if (neighbornet_objective < neighbornet_objective_minimum):
            neighbornet_objective_minimum = neighbornet_objective
            cycle_best = cycle_shifted
            
        if optimize_column_order_per_cycle and optimize_column_order:
            graphing_columns = determine_column_order(df_gather, cycle_shifted, graphing_columns, column_weights = column_weights,
                           matrix_initialization_value_column_order = matrix_initialization_value_column_order, 
                        weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, 
                        column_sorting_algorithm = column_sorting_algorithm, 
                                 verbose = verbose)
                    
    if optimize_column_order and not optimize_column_order_per_cycle:
        graphing_columns = determine_column_order(df_gather, cycle_best, graphing_columns, column_weights = column_weights,
                           matrix_initialization_value_column_order = matrix_initialization_value_column_order, 
                        weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, 
                        column_sorting_algorithm = column_sorting_algorithm, 
                                 verbose = verbose)
    return cycle_best, graphing_columns


def sort_clusters_by_agreement(df_gather, order_dict, fixed_column, reordered_column, column_weights = 'value'):
    df_temp = df_gather.copy()
    df_temp[fixed_column] = pd.Categorical(df_temp[fixed_column], ordered=True, categories = order_dict[fixed_column])
    df_temp = df_temp.sort_values(by = fixed_column).copy()
    msk = df_temp.groupby(reordered_column)[column_weights].transform('max') == df_temp[column_weights]
    reordered_list = list(df_temp.loc[msk][reordered_column].astype('str').unique())
    if len(reordered_list) < len(order_dict[reordered_column]):
        missed_clusters = [x for x in order_dict[reordered_column] if x not in reordered_list]
        order_dict[reordered_column] = reordered_list + missed_clusters
    else: 
        order_dict[reordered_column] = reordered_list 
    return order_dict


def sort_greedy_wolf(df, graphing_columns, fixed_column=None, 
                   random_initializations=1, column_weights = 'value',
                   sorting_algorithm = 'greedy_wblf', verbose = False):
    if fixed_column is None:
        fixed_column = graphing_columns[0]
    elif type(fixed_column) == int:
        fixed_column = graphing_columns[fixed_column]
        
    reordered_column = graphing_columns[graphing_columns != fixed_column]
    
    crossing_edges_objective_minimum = np.inf

    # make starting order_dict
    for i in range(random_initializations):
        order_dict = {}
        for key in graphing_columns:
            if (key == reordered_column) or (sorting_algorithm == 'greedy_wblf'):
                order_dict[key] = random.sample(list(df[key].astype('str').unique()), len(list(df[key].astype('str').unique())))
            else:
                order_dict[key] = list(df[key].astype('str').unique())
                
        temp_df = df
        if sorting_algorithm == 'greedy_wblf': 
            if i%2==0:
                order_dict = sort_clusters_by_agreement(temp_df, order_dict = order_dict,
                                           fixed_column = fixed_column,
                                          reordered_column = reordered_column, column_weights=column_weights)
                order_dict = sort_clusters_by_agreement(temp_df, order_dict = order_dict,
                                           fixed_column = reordered_column,
                                          reordered_column = fixed_column, column_weights=column_weights)
            else:
                order_dict = sort_clusters_by_agreement(temp_df, order_dict = order_dict,
                                           fixed_column = reordered_column,
                                          reordered_column = fixed_column, column_weights=column_weights)
                order_dict = sort_clusters_by_agreement(temp_df, order_dict = order_dict,
                                           fixed_column = fixed_column,
                                          reordered_column = reordered_column, column_weights=column_weights)
        elif sorting_algorithm == 'greedy_wolf':
            order_dict = sort_clusters_by_agreement(temp_df, order_dict = order_dict,
                                       fixed_column = fixed_column,
                                      reordered_column = reordered_column, column_weights=column_weights)

        if random_initializations > 1:
            crossing_edges_objective = determine_crossing_edges(df, graphing_columns, 
                                                                order_dict, 
                                                                col_weights=column_weights)
            if (crossing_edges_objective < crossing_edges_objective_minimum):
                crossing_edges_objective_minimum = crossing_edges_objective
                order_dict_best = order_dict
        else:
            order_dict_best = order_dict
    return order_dict_best

def find_colors_advanced(df_gather, graphing_columns, column_weights = 'value',
                        coloring_algorithm_advanced_option = 'leiden', resolution = 1):
    
    dist_mat, nodes = generate_distance_matrix(df_gather, graphing_columns,
                    col_weights = 'counts', matrix_initialization_value = 0, same_side_matrix_initialization_value = 0, 
                     weight_scalar = 1, nn_normalize = False)
    g = igraph.Graph.Weighted_Adjacency(dist_mat, mode = 'undirected')
    g.vs['name'] = nodes
    if coloring_algorithm_advanced_option == 'louvain':
        partition = g.community_multilevel(weights = g.es['weight'], resolution = resolution)
        mem_df = pd.DataFrame({'node' : pd.Series(nodes).str.split('~~', expand=True)[1],
                     'membership' : partition.membership})
    elif coloring_algorithm_advanced_option == 'leiden':
        partition = g.community_leiden(weights = g.es['weight'], resolution = resolution)
        mem_df = pd.DataFrame({'node' : pd.Series(nodes).str.split('~~', expand=True)[1],
                     'membership' : partition.membership})
    mem_dict = mem_df.set_index('node').to_dict()['membership']
    return mem_dict

def find_colors_reference(df_gather, graphing_columns, column_weights = 'value', reference = 'left', threshold = .5):
    if reference == 'left':
        for x in range(len(graphing_columns)-1):
            col1 = graphing_columns[x]
            col2 = graphing_columns[x+1]
                
            if x == 0:
                mem_df = pd.DataFrame({'node' : list(df_gather[col1].astype('str').unique())})
                mem_df['membership'] = range(mem_df.shape[0])
                
            df_temp = df_gather[[col1, col2, column_weights]].copy()
            df_temp = df_temp.groupby([col2, col1]).sum().reset_index()
            df_small = df_temp[df_temp.groupby([col2])[column_weights].transform(lambda x: x/sum(x)) >  threshold]
            df_small = df_small.astype('str')
            df_small = pd.merge(df_small, mem_df, left_on = col1, right_on = 'node')[[col2, 'membership']]
            df_small.columns = ['node', 'membership']
            
            for x in list(df_temp[col2].astype('str').unique()):
                if x not in df_small.node.values:
                    missing_df = pd.DataFrame({'node':[x],'membership':[np.max(mem_df.membership)+1]})
                    mem_df = pd.concat([mem_df, missing_df])
            mem_df = pd.concat([mem_df, df_small])
            mem_df['node'] = mem_df['node'].astype('str')
    elif reference == 'right':
        First = True
        for x in range(len(graphing_columns)-1)[::-1]:
            col1 = graphing_columns[x+1]
            col2 = graphing_columns[x]
                
            if First:
                mem_df = pd.DataFrame({'node' : list(df_gather[col1].astype('str').unique())})
                mem_df['membership'] = range(mem_df.shape[0])
                First = False
                
            df_temp = df_gather[[col1, col2, column_weights]].copy()
            df_temp = df_temp.groupby([col2, col1]).sum().reset_index()
            df_small = df_temp[df_temp.groupby([col2])[column_weights].transform(lambda x: x/sum(x)) >  threshold]
            df_small = df_small.astype('str')
            df_small = pd.merge(df_small, mem_df, left_on = col1, right_on = 'node')[[col2, 'membership']]
            df_small.columns = ['node', 'membership']
            
            for x in list(df_temp[col2].astype('str').unique()):
                if x not in df_small.node.values:
                    missing_df = pd.DataFrame({'node':[x],'membership':[np.max(mem_df.membership)+1]})
                    mem_df = pd.concat([mem_df, missing_df])
            mem_df = pd.concat([mem_df, df_small])
            mem_df['node'] = mem_df['node'].astype('str')
    else:  
        col1 = reference
        mem_df = pd.DataFrame({'node' : list(df_gather[col1].astype('str').unique())})
        mem_df['membership'] = range(mem_df.shape[0])
        for x in graphing_columns:
            if x != reference:
                col2 = x          
                df_temp = df_gather[[col1, col2, column_weights]].copy()
                df_temp = df_temp.groupby([col2, col1]).sum().reset_index()
                df_small = df_temp[df_temp.groupby([col2])[column_weights].transform(lambda x: x/sum(x)) >  threshold]
                df_small = df_small.astype('str')
                df_small = pd.merge(df_small, mem_df, left_on = col1, right_on = 'node')[[col2, 'membership']]
                df_small.columns = ['node', 'membership']
                
                for x in list(df_temp[col2].astype('str').unique()):
                    if x not in df_small.node.values:
                        missing_df = pd.DataFrame({'node':[x],'membership':[np.max(mem_df.membership)+1]})
                        mem_df = pd.concat([mem_df, missing_df])
                mem_df = pd.concat([mem_df, df_small])
                mem_df['node'] = mem_df['node'].astype('str')
    mem_dict = mem_df.set_index('node').to_dict()['membership']
    return mem_dict

def generate_matched_color_dict(df_gather, graphing_columns, column_weights = 'value', 
                 coloring_algorithm = 'advanced', coloring_algorithm_advanced_option = 'leiden', resolution = 1,
                cmap_name = 'tab20', threshold = .5):
    if coloring_algorithm == 'advanced':
        group_dict = find_colors_advanced(df_gather, graphing_columns, column_weights = column_weights,
                        coloring_algorithm_advanced_option = coloring_algorithm_advanced_option, resolution = resolution)
    else:
        group_dict = find_colors_reference(df_gather, graphing_columns, column_weights = column_weights, reference = coloring_algorithm,
                                          threshold = threshold)   

    cmap = plt.get_cmap(name=cmap_name, lut=len(np.unique(list(group_dict.values()))))
    available_colors = cmap(np.arange(0,cmap.N))

    color_dict = {}
    for key in group_dict.keys():
        color_dict[key] = available_colors[group_dict[key]]
    return color_dict

def plot_alluvial(df, 
                  # general function arguments
                  graphing_columns = None, column_weights = None, 
                  column1 = None, column2 = None,
                  # general sorting algorithm arguments
                sorting_algorithm = 'neighbornet', 
                  # neighbornet-specific arguments
                  optimize_column_order = True, optimize_column_order_per_cycle = False,
                  matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6,
                          weight_scalar = 5e5, 
                  # column order optimizaiton arguments
                  matrix_initialization_value_column_order = 1e6,
                  weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing",column_sorting_algorithm = "tsp", 
                  cycle_start_positions = None, 
                  # greedy wolf options
                  fixed_column = None, random_initializations = 1, 
                  # user-defined order arguments
                  order_dict=None, return_order_dict = False,
                  # alluvium arguments
                alluvium_column = None, color_alluvium = False, color_alluvium_boundary = False, alluvial_alpha = 0.5, alluvial_edge_width = 0.1,

                  # coloring algorithm arguments
                  match_colors = True,
                   coloring_algorithm = "advanced", coloring_algorithm_advanced_option = "leiden", resolution = 1, threshold = .5,
                  
                  # user-defined color arguments
                  color_dict = None, cmap_name='tab20', 
                  
                  # stratum customization options
                  color_boxes = True, include_labels_in_boxes = True, box_line_width = 1,box_width = 0.4, 
                    # stratum text_customization options
                  min_text = 4, default_text_size = 14, autofit_text = True,
                  # axis options
                 y_axis_label = False, invert_xy = False,
                    # legend options
                  include_stratum_legend = False, include_alluvium_legend = False, legend_loc = 'right', 
                    # figure size options
                  save_height = 6, save_width = 6,

                  verbose = False, savefig = False
                 ):
    df_gather = df.copy()
    figsize=(save_height, save_width)
    # Checking that necessary values are present
    if graphing_columns is None:
        if (column1 is None) or (column2 is None):
            raise ValueError("Neither graphing_columns nor both of column1, column2 are specified")
        else:
            graphing_columns = [column1, column2]
    if graphing_columns is not None:
        if len(graphing_columns) < 2:
            raise ValueError("graphing_columns must have at least 2 entries")

    if column_weights is None:
        if df.shape[1] < 2:
            raise ValueError("Dataframe must have 2 columns when column_weights is None.")
        elif df.shape[1] > 2:
            raise ValueError("column_weights must be specified when dataframe has more than 2 columns")
        else:
            df_gather = df_gather.groupby(graphing_columns).size().reset_index(name='value')
            column_weights = 'value'
    # converting values into strings
    for col in graphing_columns:
        df_gather[col] = df_gather[col].astype('str')
        
    if sorting_algorithm is not None:
        if verbose:
            print(f'Sorting Data with sorting algorithm = {sorting_algorithm}')
        if sorting_algorithm in ['greedy_wolf', 'greedy_wblf']:
            order_dict = sort_greedy_wolf(df_gather, graphing_columns, fixed_column=fixed_column, 
                   random_initializations=random_initializations, column_weights = column_weights,
                   sorting_algorithm = sorting_algorithm, verbose = verbose)
        else:
            dist_mat, nodes = generate_distance_matrix(df = df_gather, graphing_columns = graphing_columns, 
                                     col_weights=column_weights,
                      matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value,
                              weight_scalar = weight_scalar)
            if verbose:
                print(f'Sorting Distance matrix with algorithm {sorting_algorithm}')
            cycle = sort_dist_matrix(dist_mat, nodes, sorting_algorithm = sorting_algorithm)
            if verbose:
                print(f'Determining Optimal Cycle Start')
            
            cycle, graphing_columns = determine_optimal_cycle_start(df_gather, cycle, graphing_columns, column_weights=column_weights, 
                                                                    optimize_column_order = optimize_column_order,
                                     optimize_column_order_per_cycle = optimize_column_order_per_cycle, cycle_start_positions = cycle_start_positions,
                                      matrix_initialization_value_column_order = matrix_initialization_value_column_order,
                                      weight_scalar_column_order = weight_scalar_column_order, 
                                        column_sorting_metric = column_sorting_metric,column_sorting_algorithm = column_sorting_algorithm,
                                                                    verbose = verbose
                                                 )
            order_dict = get_order_dict(cycle, graphing_columns)
    if verbose:
        print("Plotting Data")

    if match_colors:
        color_dict = generate_matched_color_dict(df_gather, graphing_columns, column_weights = column_weights, 
                 coloring_algorithm = coloring_algorithm, coloring_algorithm_advanced_option = coloring_algorithm_advanced_option, 
                                  resolution = resolution, cmap_name = cmap_name)
    
    fig = plot_alluvial_internal(df_gather, graphing_columns, column_weights, 
                                 alluvium_column = alluvium_column, color_alluvium = color_alluvium, 
                                 alluvial_alpha = alluvial_alpha,color_alluvium_boundary = color_alluvium_boundary,
                                 order_dict=order_dict, color_dict = color_dict,
                                 color_boxes = color_boxes,cmap_name=cmap_name, figsize=figsize,
                                 include_labels_in_boxes = include_labels_in_boxes, box_line_width = box_line_width,box_width = box_width, 
                                 y_axis_label = y_axis_label, invert_xy = invert_xy, include_stratum_legend = include_stratum_legend, 
                                 include_alluvium_legend = include_alluvium_legend, alluvial_edge_width = alluvial_edge_width,
                                 min_text = min_text, default_text_size = default_text_size, autofit_text = autofit_text,
        )

    if type(savefig) == str:
        plt.savefig(savefig)
    
    if return_order_dict:
        return fig, order_dict
    return fig
    