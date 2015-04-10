'''
    This is the biicode hook file. It adds to your CMakeLists.txt some needed code, and copies
    another configuration file to the root folder.

    You can revert all the modifications suffered by this script in the original code. By default,
    if you're using biicode, the changes are applied:

        BII_SDSL_LITE_REVERT_CHANGES = 'False'

    To revert all, run into your current command prompt:

        Windows:    $ set BII_SDSL_LITE_REVERT_CHANGES=True
        Unix:       $ export BII_SDSL_LITE_REVERT_CHANGES=True

    and execute "bii work" to back into original form.
'''
import os
import re
import shutil
import platform


def save(path, binary_content):
    with open(path, 'wb') as handle:
        handle.write(binary_content)


def load(path):
    with open(path, 'rb') as handle:
        return handle.read()


def search_and_replace(file_path, token, replacement):
    try:
        c = load(file_path)
        c = c.replace(token, replacement)
        save(file_path, c)
    except:
        pass


def apply_changes():
    ''' Applying necessary changes to use sdsl-lite with biicode '''
    shutil.copy(os.path.join(root_folder, 'biicode', 'conf', 'biicode.conf'), root_folder)
    shutil.copy(os.path.join(root_folder, 'biicode', 'conf', 'ignore.bii'), root_folder)
    search_and_replace(cmakelist_path, cmakelist_token, cmakelist_replacement)


def revert_changes():
    ''' Revert all the biicode changes code '''
    os.remove(biicode_conf_path)
    os.remove(ignore_path)
    if os_platform == "Windows":
        search_and_replace(cmakelist_path, cmakelist_replacement_win, cmakelist_token)
    else:
        search_and_replace(cmakelist_path, cmakelist_replacement, cmakelist_token)


# Main code
os_platform = platform.system()
BII_SDSL_LITE_REVERT_CHANGES = os.environ.get('BII_SDSL_LITE_REVERT_CHANGES', 'False')

root_folder = bii.block_folder if os.path.exists(bii.block_folder) else bii.project_folder
biicode_conf_path = os.path.join(root_folder, 'biicode.conf')
ignore_path = os.path.join(root_folder, 'ignore.bii')

cmakelist_path = os.path.join(root_folder, "CMakeLists.txt")
cmakelist_token = "#COMMENT REPLACED BY BIICODE"
cmakelist_replacement = '''if(BIICODE)
include(biicode/cmake/biicode.cmake)
return()
endif()'''
cmakelist_replacement_win = "if(BIICODE)\r\ninclude(biicode/cmake/biicode.cmake)\r\nreturn()\r\nendif()"

try:
    # Apply or revert changes
    if BII_SDSL_LITE_REVERT_CHANGES == 'False':
        if "if(BIICODE)" in load(cmakelist_path):
            print "Hook: changes just applied"
        else:
            print "Hook: applying changes"
            apply_changes()
    else:
        if "if(BIICODE)" not in load(cmakelist_path):
           print "Hook: changes just reverted"
        else:
            print "Hook: reverting changes"
            revert_changes()
except Exception as e:
    print "Exception: %s" % e

