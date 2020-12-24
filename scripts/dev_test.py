import time
from glowingmeme.build_data.build_dataset import BuildDatasetCVA, BuildDatasetCipapi

start_time = time.time()
bd_cva = BuildDatasetCVA()

bd_cipapi = BuildDatasetCipapi(bd_cva.build_dataset())
bd_cipapi.build_dataset()
bd_cipapi.save_data_to_csv("hello_there.csv")

print("--- %s seconds ---" % (time.time() - start_time))
