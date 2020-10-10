from glowingmeme.build_data.build_dataset import BuildDatasetCVA, BuildDatasetCipapi

bd_cva = BuildDatasetCVA()
main_dataset = bd_cva.build_dataset()

bd_cipapi = BuildDatasetCipapi(main_dataset)
main_dataset = bd_cipapi.build_dataset()

x = 1
