from glowingmeme.build_data.build_dataset import BuildDatasetCVA

bd = BuildDatasetCVA()

x, y = bd._query_cva_archived_positive_cases()
h = 1