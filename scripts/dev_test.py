import time
from glowingmeme.build_data.build_dataset_cva import BuildDatasetCVA
from glowingmeme.build_data.build_dataset_cipapi import  BuildDatasetCipapi
from glowingmeme.build_data.build_dataset_cellbase import BuildDatasetCellbase

start_time = time.time()
bd_cva = BuildDatasetCVA()
bd_cva.build_dataset()

bd_cipapi = BuildDatasetCipapi(bd_cva.main_dataset)
bd_cipapi.build_dataset()

bd_cellbase = BuildDatasetCellbase(bd_cipapi.main_dataset)
bd_cellbase.build_dataset()

bd_cellbase.save_data_to_csv("hello_there.csv")

print("--- %s seconds ---" % (time.time() - start_time))
