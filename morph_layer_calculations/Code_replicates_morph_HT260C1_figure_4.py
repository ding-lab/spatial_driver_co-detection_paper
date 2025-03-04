# Code related to running Morph for figure 4
# HT260C1-Th1K1L1U1
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT260C1_rerun
# conda activate spatial_driver_layer
import csv
import gzip
import matplotlib.pyplot
import numpy
import skimage
x_location = []
y_location = []
feature_name = []
cell_id = []

filename = '/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT260C1-Th1K1L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz'
annotation = '/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v2/HT260C1-Th1K1L1U1/Manual_individual_cell_type_v3_subclone.tsv'

groups = {}
with open(annotation, 'r') as f:
    dict_reader = csv.DictReader(f, delimiter='\t')
    for line in dict_reader:
        groups[line['cell_id']] = line['Manual_individual_cell_type_v3_sub']
    
with gzip.open(filename, 'rt') as f:
    reader = csv.DictReader(f)
    for row in reader:
        x_location.append(float(row['x_location']))
        y_location.append(float(row['y_location']))
        feature_name.append(row['feature_name'])
        cell_id.append(row['cell_id'])

density = 10
testX = {}
tempSet = []

for c in set(groups.values()):
    testX[c] = numpy.zeros((int(max(x_location) / density + 1), int(max(y_location) / density + 1)))
for i in range(len(x_location)):
    if cell_id[i] in groups.keys() and groups[cell_id[i]]:
        testX[groups[cell_id[i]]][int(float(x_location[i]) / density), int(float(y_location[i]) / density)] += 1

centroid = '/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT260C1-Th1K1L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz'

x_centroid = {}
y_centroid = {}
with gzip.open(centroid, 'rt') as f:
    reader = csv.DictReader(f)
    for row in reader:
        x_centroid[row['cell_id']] = row['x_centroid']
        y_centroid[row['cell_id']] = row['y_centroid']

testX2 = {}

for c in set(groups.values()):
    print(c)
    testX2[c] = numpy.zeros((int(max(x_location) / density + 1), int(max(y_location) / density + 1)))
for i in x_centroid.keys():
    if i in groups.keys():
        testX2[groups[i]][int(float(x_centroid[i]) / density), int(float(y_centroid[i]) / density)] += 1

def layering(image, element=None, method='min', tissue=1):
    image += 1 - tissue
    eroded = skimage.morphology.erosion(image, element)
    layers = numpy.multiply(image - eroded, tissue)
    while not (numpy.multiply(image, tissue) == numpy.multiply(eroded, tissue)).all():
        image = eroded
        eroded = skimage.morphology.erosion(image, element)
        layers = layers + numpy.multiply(image - eroded, tissue) * (numpy.max(layers) + 1)
    return layers

def hexagon(side):
    return numpy.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]])

tumor_1 = testX['Tumor_1'] + 0
tumor_1 += testX2['Tumor_1']
tumor_2 = testX['Tumor_2'] + 0
tumor_2 += testX2['Tumor_2']

# [100:440, 290:]
tumor_1_filtered = numpy.zeros_like(testX['Tumor_1'])
tumor_1_filtered[290:, 100:440] = tumor_1[290:, 100:440]
tumor_2_filtered = tumor_2.copy()
tumor_2_filtered[290:, 100:440] = 0

mask = numpy.maximum.reduce([testX[c] for c in testX.keys()]).T > 0
temp1 = layering(1 - (tumor_1_filtered > 0).astype(int), tissue=numpy.ones_like(tumor_1_filtered).astype(int))
temp1 = temp1.T
temp1 = numpy.multiply(temp1, mask)


temp2 = layering(1 - (tumor_2_filtered > 0).astype(int), tissue=numpy.ones_like(tumor_2_filtered).astype(int))
temp2 = temp2.T
temp2 = numpy.multiply(temp2, mask)

matplotlib.pyplot.figure(figsize=(20, 20))
matplotlib.pyplot.subplot(1, 2, 1)
matplotlib.pyplot.imshow(temp1)
matplotlib.pyplot.colorbar()
matplotlib.pyplot.subplot(1, 2, 2)
matplotlib.pyplot.imshow(temp2)
matplotlib.pyplot.colorbar()
matplotlib.pyplot.savefig('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT260C1_rerun/HT260_layer_morph_contour.pdf')
matplotlib.pyplot.clf()

cells = {}
for c in set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_1']):
    bins = []
    for layer in range(20):
        bins.append(sum(testX2[c].T[(3 * layer < temp2) & (temp2 <= 3 * (layer + 1))].flatten()))
    cells[c] = bins

with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT260C1_rerun/HT260_layer_stacked_bar_plot_tumor_2_source.csv', 'w') as f:
    dict_writer = csv.DictWriter(f, fieldnames=['cell_type', 'count'])
    dict_writer.writeheader()
    for c in cells.keys():
        dict_writer.writerow({'cell_type': c, 'count': cells[c]})

bins = [str(i) + '-' + str(j) for i, j in zip(range(0, 30 * 20, 30), range(30, 30 * 21, 30))]
header = ['cell_type']
[header.append(i) for i in bins]

field_dict = {}
for i,v in enumerate(header):
    field_dict[v] = i

for i in range(len(cells['Macrophage'])):
    proportion = 0
    for c in set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_1', 'Macrophage_LowCount', 'Fibroblast_LowCount']):
        proportion += cells[c][i]
    for c in set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_1', 'Macrophage_LowCount', 'Fibroblast_LowCount']):
        cells[c][i] = cells[c][i] / proportion

matplotlib.pyplot.figure(figsize=(6, 5))
bottom = numpy.zeros(20)
for c, color in zip(set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_1', 'Macrophage_LowCount', 'Fibroblast_LowCount']), matplotlib.colormaps['tab20'].colors):
    matplotlib.pyplot.bar([str(i) + ' - ' + str(j) for i, j in zip(range(0, 30 * 20, 30), range(30, 30 * 21, 30))], cells[c], label=c, bottom=bottom, color=color)
    bottom += cells[c]
matplotlib.pyplot.gca().legend(bbox_to_anchor=(1.2, 1.0))
_ = matplotlib.pyplot.xticks(rotation=90)
matplotlib.pyplot.savefig('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT260C1_rerun/HT260_layer_stacked_bar_plot_tumor_2_source.pdf')
matplotlib.pyplot.clf()

cells = {}
for c in set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_2']):
    bins = []
    for layer in range(20):
        bins.append(sum(testX2[c].T[(3 * layer < temp1) & (temp1 <= 3 * (layer + 1))].flatten()))
    cells[c] = bins

with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT260C1_rerun/HT260_layer_stacked_bar_plot_tumor_1_source.csv', 'w') as f:
    dict_writer = csv.DictWriter(f, fieldnames=['cell_type', 'count'])
    dict_writer.writeheader()
    for c in cells.keys():
        dict_writer.writerow({'cell_type': c, 'count': cells[c]})

with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT260C1_rerun/HT260_layer_stacked_bar_plot_tumor_2_source.tsv', 'w') as f:
    dict_writer = csv.DictWriter(f, fieldnames=header, delimiter='\t')
    dict_writer.writeheader()
    for c in cells.keys():
        newline_dict = {}
        for field_key in field_dict.keys():
            field_num = field_dict[field_key]
            if field_num == 0:
                newline_dict[field_key] = c
            else:
                newline_dict[field_key] = cells[c][field_num - 1]
        dict_writer.writerow(newline_dict)

for i in range(len(cells['Macrophage'])):
    proportion = 0
    for c in set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_2', 'Macrophage_LowCount', 'Fibroblast_LowCount']):
        proportion += cells[c][i]
    for c in set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_2', 'Macrophage_LowCount', 'Fibroblast_LowCount']):
        cells[c][i] = cells[c][i] / proportion

matplotlib.pyplot.figure(figsize=(6, 5))
bottom = numpy.zeros(20)
for c, color in zip(set(testX.keys()) - set(['Necrosis?', 'LowCount', 'Tumor_2', 'Macrophage_LowCount', 'Fibroblast_LowCount']), matplotlib.colormaps['tab20'].colors):
    matplotlib.pyplot.bar([str(i) + ' - ' + str(j) for i, j in zip(range(0, 30 * 20, 30), range(30, 30 * 21, 30))], cells[c], label=c, bottom=bottom, color=color)
    bottom += cells[c]
matplotlib.pyplot.gca().legend(bbox_to_anchor=(1.2, 1.0))
_ = matplotlib.pyplot.xticks(rotation=90)
matplotlib.pyplot.savefig('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT260C1_rerun/HT260_layer_stacked_bar_plot_tumor_1_source.pdf')
