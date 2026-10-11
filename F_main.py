
from F1_static import F1_static
from FD31_DynamicCore import FD31_DynamicCore
from FD32_postprocessing import FD32_postprocessing


# ============================================================================
# USAGE
# ============================================================================
if __name__ == "__main__":

    for i in range (0,5,1):

        FileToOpen = f'primer_{i}.txt'  # Should be file name+extension, to read and write new results
        print(f"FileToOpen", FileToOpen)

        F1_static(FileToOpen)

        result = FD31_DynamicCore(FileToOpen)  # calculate independent coordinates, and store them into more file intervals
        print(f"Completed {result} intervals")

        FD32_postprocessing(FileToOpen)
        print(f"Completed Postprocessing")
    print('End of all')





