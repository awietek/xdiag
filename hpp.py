#!/usr/bin/env python
import os

root = "."

# # rename .h -> .hpp
# for path, subdirs, files in os.walk(root):
#     for name in [name for name in files if name[-2:] == ".h"]:
#         orig = os.path.join(path, name)
#         dest = os.path.join(path, name + "pp")
#         print("renaming")
#         print(orig)
#         print(dest)

#         input("Press Enter to continue...")
#         os.rename(orig, dest)
#         print()
        

# fix the includes
for path, subdirs, files in os.walk(root):
    for name in [name for name in files if name[-4:] == ".hpp" or name[-4:] == ".cpp"]:
        orig = os.path.join(path, name)
        dest = os.path.join(path, name + "pp")
        print("At file")
        print(orig)

        with open(orig, 'r') as fl:
            for line in fl.readlines():
                if "#include" in line and (".h\"" in line or ".h>" in line):
                    print(line)

        
        input("Press Enter to continue...")

        print("-------------------")
        print()
