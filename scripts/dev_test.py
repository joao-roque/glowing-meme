import time
from guppy import hpy
from glowingmeme.build_data.build_dataset import BuildDatasetCVA, BuildDatasetCipapi

memory_tracking = hpy()

start_time = time.time()
bd_cva = BuildDatasetCVA()
main_dataset = bd_cva.build_dataset()

bd_cipapi = BuildDatasetCipapi(main_dataset)
bd_cipapi.build_dataset()
bd_cipapi.save_data_to_csv("hello_there.csv")

print("--- %s seconds ---" % (time.time() - start_time))
print(memory_tracking.heap())

x = 1
