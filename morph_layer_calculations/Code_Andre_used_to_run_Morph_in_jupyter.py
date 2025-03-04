# Code used to run Morph for layer calculations on each sample then used as input for figures 5 and 6
# stuff related to running just a single sample has been commented out.

import csv
import gzip
import matplotlib.pyplot
import numpy

filenames = [
'/diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT227P1-S1H1L1U1__20240411__214815/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT242P1-S1H4L4U1__20240411__214815/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__SP001P1-Fp1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__HT270P1-S1H1A1US2_1__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__SP002C1-Fp1U2__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT060P1-S1R1Fp1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L4U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT125P1-S1H4A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT125P1-S1H8A1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT179C1-T1Fp3L5U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT185P1-S1H2L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/transcripts.csv.gz',
]

annotations = [
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT227P1-S1H1L1U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT242P1-S1H4L4U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v1/SP001P1-Fp1U1/Manual_individual_cell_type_v1.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT270P1-S1H1A1US2_1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v1/SP002C1-Fp1U2/Manual_individual_cell_type_v1.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT060P1-S1R1Fp1U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT061P1-S1P1A1L1U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT061P1-S1P1A1L4U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT125P1-S1H4A1L1U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT125P1-S1H8A1U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v1/HT179C1-T1Fp3L5U1/HT179C1-T1Fp3L5U1_cell_type_col_fix.tsv',
'/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT185P1-S1H2L1U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv',
]

centroids = [
'/diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT227P1-S1H1L1U1__20240411__214815/cells.csv.gz',
'/diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT242P1-S1H4L4U1__20240411__214815/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__SP001P1-Fp1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__HT270P1-S1H1A1US2_1__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__SP002C1-Fp1U2__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT060P1-S1R1Fp1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L4U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT125P1-S1H4A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT125P1-S1H8A1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT179C1-T1Fp3L5U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
'/diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT185P1-S1H2L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/cells.csv.gz',
]

heads = [
'group',
'group',
'Manual_individual_cell_types_v1',
'group',
'Manual_individual_cell_types_v1',
'group',
'group',
'group',
'group',
'group',
'Manual_individual_cell_types_v1',
'group',
]


x_centroids = {}
y_centroids = {}
cell_ids = {}
testX = {}
testX2 = {}
#for filename, annotation, centroid, head in zip([filenames[0]], [annotations[0]], [centroids[0]], [heads[0]]):
for filename, annotation, centroid, head in zip(filenames, annotations, centroids, heads):
    print(annotation)
    x_location = []
    y_location = []
    feature_name = []
    cell_id = []
    groups = {}
    with open(annotation, 'r') as f:
        dict_reader = csv.DictReader(f, delimiter='\t')
        for line in dict_reader:
            groups[line['cell_id']] = line[head]
    with gzip.open(filename, 'rt') as f:
        reader = csv.DictReader(f)
        for row in reader:
            x_location.append(float(row['x_location']))
            y_location.append(float(row['y_location']))
            feature_name.append(row['feature_name'])
            cell_id.append(row['cell_id'])
    density = 10
    testX[annotation] = {}
    for c in set(groups.values()):
        testX[annotation][c] = numpy.zeros((int(max(x_location) / density + 1), int(max(y_location) / density + 1)))
    for i in range(len(x_location)):
        if cell_id[i] in groups.keys() and groups[cell_id[i]]:
            testX[annotation][groups[cell_id[i]]][int(float(x_location[i]) / density), int(float(y_location[i]) / density)] += 1
    x_centroids[annotation] = {}
    y_centroids[annotation] = {}
    cell_ids[annotation] = []
    with gzip.open(centroid, 'rt') as f:
        reader = csv.DictReader(f)
        for row in reader:
            x_centroids[annotation][row['cell_id']] = row['x_centroid']
            y_centroids[annotation][row['cell_id']] = row['y_centroid']
            cell_ids[annotation].append(row['cell_id'])
    testX2[annotation] = {}
    for c in set(groups.values()):
        print(c)
        testX2[annotation][c] = numpy.zeros((int(max(x_location) / density + 1), int(max(y_location) / density + 1)))
    for i in x_centroids.keys():
        if i in groups.keys():
            testX2[annotation][groups[i]][int(float(x_centroid[i]) / density), int(float(y_centroid[i]) / density)] += 1

import skimage
import pandas
# this function is a re-implementation of the morph layer calculation on Xenium data: https://github.com/ding-lab/morph
def layering(image, element=None, method='min', tissue=1):
    image += 1 - tissue
    eroded = skimage.morphology.erosion(image, element)
    layers = numpy.multiply(image - eroded, tissue)
    while not (numpy.multiply(image, tissue) == numpy.multiply(eroded, tissue)).all():
        image = eroded
        eroded = skimage.morphology.erosion(image, element)
        layers = layers + numpy.multiply(image - eroded, tissue) * (numpy.max(layers) + 1)
    return layers

sample = 'HT227P1-S1H1L1U1'

input_table = pandas.read_csv("/diskmnt/Projects/Users/andretargino/Draft/andre/Morph_input_table_v7.txt", sep='\t', header=0)
sample_list = list(input_table["Sample_ID"])
print(sample_list)
#for sample_id in sample_list:
###sample_list = list(testX.keys())

sample_ids = [
'HT227P1-S1H1L1U1',
'HT242P1-S1H4L4U1',
'SP001P1-Fp1U1',
'HT270P1-S1H1A1US2_1',
'SP002C1-Fp1U2',
'HT060P1-S1R1Fp1U1',
'HT061P1-S1P1A1L1U1',
'HT061P1-S1P1A1L4U1',
'HT125P1-S1H4A1L1U1',
'HT125P1-S1H8A1U1',
'HT179C1-T1Fp3L5U1',
'HT185P1-S1H2L1U1',
'HT403C1-S1H1A1U1',
]

sample_list.remove("HT403C1-S1H1A1U1")

#for sample_id, annotation in zip(['HT061P1-S1P1A1L1U1'], ['/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT061P1-S1P1A1L1U1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv']):
#for sample_id, annotation in zip(['HT270P1-S1H1A1US2_1'], ['/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/HT270P1-S1H1A1US2_1_PDAC_merge_primary_KRAS_20241205_single_SCT_cell_type_merge_v3.tsv']):
#for sample_id, sample in zip(['HT227P1-S1H1L1U1'], ['output-XETG00122__0022792__HT227P1-S1H1L1U1__20240411__214815']):#['HT270P1-S1H1A1US2_1']:#sample_list:
for sample_id, annotation in zip(sample_ids, annotations):
    sample = annotation
#    sample_id = sample
    tumor_csv = input_table.loc[input_table.index[input_table.loc[:,'Sample_ID'] == sample_id],["Xenium_neoplastic_labels_v7"]].values[0][0]
    tumor_label_list = tumor_csv.strip().split(',')
    if 'PanIN' in tumor_label_list:
        tumor_2_label_list = ['PanIN']
        tumor_1_label_list = tumor_label_list.copy()
        tumor_1_label_list.remove('PanIN')
    else:
        tumor_2_label_list = []
        tumor_1_label_list = tumor_label_list.copy()
    unknown_csv = input_table.loc[input_table.index[input_table.loc[:,'Sample_ID'] == sample_id],["Xenium_unknown_low_quality_labels_v7"]].values[0][0]
#    print(unknown_csv)
    if unknown_csv in [numpy.NaN, 'nan']:
        unknown_label_list = []
    else:
        unknown_label_list = unknown_csv.strip().split(',')
    print(sample)
    print(tumor_1_label_list)
    print(tumor_2_label_list)
    print(unknown_label_list)
    tumor_1 = testX[sample]['Tumor'] + 0
    tumor_1 += testX2[sample]['Tumor']
    if len(tumor_1_label_list) > 1:
        tumor_1 += testX[sample][tumor_1_label_list[-1]]
        tumor_1 += testX2[sample][tumor_1_label_list[-1]]
    tumor_2 = numpy.zeros_like(testX[sample]['Tumor'])
    if len(tumor_2_label_list):
        tumor_2 = testX[sample]['PanIN'] + 0
        tumor_2 += testX2[sample]['PanIN']
#    tumor_1 = numpy.zeros_like(testX[sample]['Tumor'])
    if sample_id == 'HT227P1-S1H1L1U1':
        tumor_1 = testX[sample]['Tumor'] + 0
        tumor_1 += testX2[sample]['Tumor']
        tumor_1[340:380, 210:230] = 0
        tumor_1[360:380, 193:230] = 0
        tumor_2 = numpy.zeros_like(testX[sample]['PanIN'])
        tumor_2[340:380, 210:230] = testX[sample]['PanIN'][340:380, 210:230]
        tumor_2[340:380, 210:230] += testX2[sample]['PanIN'][340:380, 210:230]
        tumor_2[360:380, 193:230] = testX[sample]['PanIN'][360:380, 193:230]
        tumor_2[360:380, 193:230] += testX2[sample]['PanIN'][360:380, 193:230]
    if sample_id == 'HT270P1-S1H1A1US2_1':
        tumor_1 = testX[sample]['Tumor'] + 0
        tumor_1 += testX2[sample]['Tumor']
        tumor_1[165:210, 179:220] = 0
        tumor_2 = testX[sample]['PanIN'] + 0
        tumor_2 += testX2[sample]['PanIN']
        tumor_2[165:210, 179:220] = 0
    if sample_id == 'HT061P1-S1P1A1L1U1':
        tumor_2 = testX[sample]['PanIN'] + 0
        tumor_2 += testX2[sample]['PanIN']
        tumor_2[70:120, 330:370] = 0
        tumor_2[497:555, 134:292] = 0
#        tumor_2 = numpy.zeros_like(testX[sample]['PanIN'])
#        tumor_2[70:120, 330:370] = testX[sample]['PanIN'][70:120, 330:370] 
#        tumor_2[70:120, 330:370] += testX2[sample]['PanIN'][70:120, 330:370]
#        tumor_2[497:555, 134:292] = testX[sample]['PanIN'][497:555, 134:292]
#        tumor_2[497:555, 134:292] += testX2[sample]['PanIN'][497:555, 134:292]
    matplotlib.pyplot.figure(figsize=(20, 20))
    matplotlib.pyplot.title(sample)
    matplotlib.pyplot.subplot(1, 2, 1)
    matplotlib.pyplot.imshow(tumor_1)
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.subplot(1, 2, 2)
    matplotlib.pyplot.imshow(tumor_2)
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.savefig('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/rerun_check_morph_andre_v7/' + sample_id + '_tumor_layer_origin.pdf')
    mask = numpy.maximum.reduce([testX[sample][c] for c in testX[sample].keys()]) > 0
    temp1 = layering(1 - (tumor_1 > 0).astype(int), tissue=numpy.ones_like(tumor_1).astype(int))
    temp1 = numpy.multiply(temp1, mask)
    temp2 = layering(1 - (tumor_2 > 0).astype(int), tissue=numpy.ones_like(tumor_2).astype(int))
    temp2 = numpy.multiply(temp2, mask)
    matplotlib.pyplot.figure(figsize=(20, 20))
    matplotlib.pyplot.title(sample)
    matplotlib.pyplot.subplot(1, 2, 1)
    matplotlib.pyplot.imshow(temp1)
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.subplot(1, 2, 2)
    matplotlib.pyplot.imshow(temp2)
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.savefig('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/rerun_check_morph_andre_v7/' + sample_id + '_layer_morph_contour.pdf')
    cells = {}
    for c in set(testX[sample].keys()) - set(unknown_label_list):
        bins = []
        for layer in range(numpy.max(temp1)):
            bins.append(sum(testX2[sample][c][temp1 == layer].flatten()))
        cells[c] = bins
    with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/rerun_check_morph_andre_v7/' + sample_id + '_layer_stacked_bar_plot_tumor_1_source.csv', 'w') as f:
        dict_writer = csv.DictWriter(f, fieldnames=['cell_type'] + [str(i) for i in range(numpy.max(temp1))])
        dict_writer.writeheader()
        for c in cells.keys():
            d = {}
            d['cell_type'] = c
            for i in range(numpy.max(temp1)):
                d[str(i)] = int(cells[c][i])
            dict_writer.writerow(d)
    with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/rerun_check_morph_andre_v7/' + sample_id + '_cell_layer_tumor_1_source.csv', 'w') as f:
        dict_writer = csv.DictWriter(f, fieldnames=['cell_id', 'layer'])
        dict_writer.writeheader()
        for c in cell_ids[sample]:
            dict_writer.writerow({'cell_id': c, 'layer': temp1[int(float(x_centroids[sample][c]) / density), int(float(y_centroids[sample][c]) / density)]})
    cells = {}
    for c in set(testX[sample].keys()) - set(unknown_label_list):
        bins = []
        for layer in range(numpy.max(temp2)):
            bins.append(sum(testX2[sample][c][temp2 == layer].flatten()))
        cells[c] = bins
    with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/rerun_check_morph_andre_v7/' + sample_id + '_layer_stacked_bar_plot_tumor_2_source.csv', 'w') as f:
        dict_writer = csv.DictWriter(f, fieldnames=['cell_type'] + [str(i) for i in range(numpy.max(temp2))])
        dict_writer.writeheader()
        for c in cells.keys():
            d = {}
            d['cell_type'] = c
            for i in range(numpy.max(temp2)):
                d[str(i)] = int(cells[c][i])
            dict_writer.writerow(d)
    with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/rerun_check_morph_andre_v7/' + sample_id + '_cell_layer_tumor_2_source.csv', 'w') as f:
        dict_writer = csv.DictWriter(f, fieldnames=['cell_id', 'layer'])
        dict_writer.writeheader()
        for c in cell_ids[sample]:
            dict_writer.writerow({'cell_id': c, 'layer': temp2[int(float(x_centroids[sample][c]) / density), int(float(y_centroids[sample][c]) / density)]})


# The code used for the tumor microregion calculation and 10µm expansion in Extended Data Figure5c,d for HT284P1 calculation
import csv
import gzip
import matplotlib.pyplot
import numpy
import scipy
import skimage
import tifffile
import subprocess

cells = '/diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__HT284P1-S1H1A1U1__20240417__211045/cells.csv.gz'
h_and_e = '/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT284P1-S1H1A1U1/v9_bs0m60/he_aligned.ome.tif'
tumor = '/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/clinical_annotations/HT284P1-S1H1A1U1/v9_s0m60/HT284P1-S1H1A1U1_v9_bs0m60_he_aligned_annotation_updated_tumor_only_20240918.ome.tif'
necrosis = '/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/clinical_annotations/HT284P1-S1H1A1U1/v9_s0m60/HT284P1-S1H1A1U1_v9_bs0m60_he_aligned_annotation_updated_necrosis_only_20240918.ome.tif'
stroma = '/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/clinical_annotations/HT284P1-S1H1A1U1/v9_s0m60/HT284P1-S1H1A1U1_v9_bs0m60_he_aligned_annotation_updated_stroma_recovery_20240918.ome.tif'

cell_id = []
x_centroid = []
y_centroid = []
with gzip.open(cells, 'rt') as f:
    reader = csv.DictReader(f)
    for row in reader:
        cell_id.append(row['cell_id'])
        x_centroid.append(float(row['x_centroid']))
        y_centroid.append(float(row['y_centroid']))

with tifffile.TiffFile(tumor) as f:
    tumor = f.series[0].levels[0].asarray()

with tifffile.TiffFile(necrosis) as f:
    necrosis = f.series[0].levels[0].asarray()

with tifffile.TiffFile(stroma) as f:
    stroma = f.series[0].levels[0].asarray()

gray = skimage.color.rgb2gray(necrosis)
binary = gray > 50 / 255
filled = scipy.ndimage.binary_fill_holes(binary)
filtered = skimage.morphology.remove_small_objects(filled, 50)
distance = scipy.ndimage.distance_transform_edt(filtered)
# eroded = skimage.morphology.erosion(distance,skimage.morphology.square(3))
# eroded = numpy.zeros_like(distance)
eroded = skimage.morphology.erosion(distance,skimage.morphology.square(1))
eroded[distance > 10 / 0.2125] = 1
label, _ = scipy.ndimage.label(eroded, numpy.ones((3, 3)))
with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT284P1_microregion_10um_expansion_rerun/testHT284P1necrosisv2.csv', 'w') as f:
    dict_writer = csv.DictWriter(f, fieldnames=['barcode', 'label'])
    dict_writer.writeheader()
    for i in range(len(cell_id)):
        if label[int(y_centroid[i] / 0.2125), int(x_centroid[i] / 0.2125)]:
            dict_writer.writerow({'barcode': cell_id[i], 'label': label[int(y_centroid[i] / 0.2125), int(x_centroid[i] / 0.2125)]})

gray = skimage.color.rgb2gray(stroma)
# gray = skimage.color.rgb2gray(tumor)
binary = gray > 50 / 255
filled = scipy.ndimage.binary_fill_holes(binary)
filtered = skimage.morphology.remove_small_objects(filled, 50)
distance = scipy.ndimage.distance_transform_edt(filtered)
eroded = numpy.zeros_like(distance)
# we were aiming to include cells within 10 microns of tumor 10/0.2125 ~ 47.05
# eroded = skimage.morphology.dilation(distance,skimage.morphology.disk(47))
eroded[distance > 10 / 0.2125] = 1
label, _ = scipy.ndimage.label(eroded, numpy.ones((3, 3)))
with open('/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/HT284P1_microregion_10um_expansion_rerun/testHT284P1stromav2.csv', 'w') as f:
    dict_writer = csv.DictWriter(f, fieldnames=['barcode', 'label'])
    dict_writer.writeheader()
    for i in range(len(cell_id)):
        if label[int(y_centroid[i] / 0.2125), int(x_centroid[i] / 0.2125)]:
            dict_writer.writerow({'barcode': cell_id[i], 'label': label[int(y_centroid[i] / 0.2125), int(x_centroid[i] / 0.2125)]})
