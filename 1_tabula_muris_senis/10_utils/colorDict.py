import matplotlib.pyplot as plt

def tissue_color_dict():
    tissues = ['Bladder', 'Brain_Myeloid', 'Brain_Non-Myeloid', 'Fat', 'Heart', 'Kidney', 'Large_Intestine', 'Limb_Muscle', 'Liver', 'Lung', 'Mammary_Gland', 'Marrow', 'Pancreas', 'Skin', 'Spleen', 'Thymus', 'Tongue', 'Trachea', 'Aorta', 'Diaphragm', 'BAT', 'GAT', 'MAT', 'SCAT']
    tissues.sort()
    
    import matplotlib.colors as pltcolors
    
    cmap = [plt.cm.get_cmap("tab20b")(0)] # Aorta
    for i in range(3,5): # BAT, Bladder
        cmap.append(plt.cm.get_cmap("tab20b")(i))
    for i in range(6,9): # Brain_Myeloid, Brain_Non_Myeloid, Diaphgram
        cmap.append(plt.cm.get_cmap("tab20b")(i))
    for i in range(17,18): # Fat
        cmap.append(plt.cm.get_cmap("tab20c")(i))
    for i in range(9,13): # GAT, Heart, Kidney, Large_Intestine
        cmap.append(plt.cm.get_cmap("tab20b")(i))
    for i in range(14,20): # Limb_Muscle, Liver, Lung, MAT, Mammary_Gland, Marrow
        cmap.append(plt.cm.get_cmap("tab20b")(i))
    for i in range(0,2): # Pancreas, SCAT
        cmap.append(plt.cm.get_cmap("tab20c")(i))
    for i in range(4,6): # Skin, Spleen
        cmap.append(plt.cm.get_cmap("tab20c")(i))
    for i in range(8,9): # Thymus
        cmap.append(plt.cm.get_cmap("tab20c")(i))
    for i in range(10,11): # Tongue
        cmap.append(plt.cm.get_cmap("tab20c")(i))
    for i in range(16,17): # Trachea 
        cmap.append(plt.cm.get_cmap("tab20c")(i))
    
    color_dict = {}
    j=0
    for t in tissues:
        color_dict[t] = pltcolors.to_hex(cmap[j])
        j+=1
    
    return color_dict

def plot_colortable(colors, title, sort_colors=True, emptycols=0):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    topmargin = 40

    # Sort colors by hue, saturation, value and name.
    by_hsv = [(v, k) for k, v in colors.items()]
    
    if sort_colors is True:
        by_hsv = sorted(by_hsv)
    names = [name for hsv, name in by_hsv]

    n = len(names)
    ncols = 4 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24, loc="left", pad=10)

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=colors[name], linewidth=18)

    return fig

def age_color_dict():
    ages = ['1m','3m','18m','21m','24m','30m']
    
    import matplotlib.colors as pltcolors
    
    cmap = plt.cm.get_cmap("YlGnBu")
        
    age_color_dict = {}
    j=1/len(ages)
    for a in ages:
        age_color_dict[a] = pltcolors.to_hex(cmap(j))
        j+=1/len(ages)
        

    
    return age_color_dict

def method_color_dict():
    methods = ['droplet','facs']
    
    import matplotlib.colors as pltcolors
    
    method_color_dict = {}
    method_color_dict['droplet'] = '#90ee90'
    method_color_dict['facs'] = '#006400'
    
    return method_color_dict

