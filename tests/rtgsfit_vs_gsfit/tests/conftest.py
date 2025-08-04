import os
import subprocess

import pytest
import MDSplus
import mdsthin


from rtgsfit_vs_gsfit import cnst, initialise_rtgsfit_node, replay_gsfit

@pytest.fixture(scope="session")
def run_gsfit_and_save_dict():
    """
    Run GSFIT, then save the results to a dictionary.
    """
    replay_gsfit.replay_gsfit()

@pytest.fixture(scope="session")
def rtgsfit_mds_nodegen():
    """
    Delete the existing RTGSFIT MDSplus node and ensure it is removed.
    """

    # Delete the existing RTGSFIT MDSplus node
    def delete_node_recursive(tree, node):
        for child in node.getNodeWild('*'):
            delete_node_recursive(tree, child)
        print(f"Deleting node {node.getPath()}")
        tree.deleteNode(node.getPath())
    tree = MDSplus.Tree("RTGSFIT", cnst.PULSE_NUM_WRITE, "EDIT")
    try:
        node = tree.getNode(f":{cnst.RUN_NAME}")
        delete_node_recursive(tree, node)
        print(f"Deleted node {cnst.RUN_NAME} and its children.")
    except Exception as e:
        print(f"Failed to delete node {cnst.RUN_NAME}: {e}")
    tree.write()
    tree.close()

    # Re-open tree to check if node exists (usually best to re-open fresh)
    def node_exists(tree, path):
        try:
            tree.getNode(path)
            return True
        except Exception:
            return False
    tree_check = MDSplus.Tree("RTGSFIT", cnst.PULSE_NUM_WRITE)
    if not node_exists(tree_check, f":{cnst.RUN_NAME}"):
        deleted = True
    else:
        deleted = False
    assert deleted, f"Node {cnst.RUN_NAME} was not deleted successfully."

    # # Initialise the RTGSFIT MDSplus node
    initialise_rtgsfit_node.initialise_rtgsfit_node()

@pytest.fixture(scope="session")
def compile_rtgsfit(rtgsfit_mds_nodegen):
    """
    Clean and then compile RTGSFIT.
    """

    # Clean the RTGSFIT source directory
    os.chdir(cnst.RTGSFIT_SRC_PATH)
    subprocess.run(["make", "clean"], cwd=cnst.RTGSFIT_SRC_PATH, check=True)
    # Check if the object files are removed
    object_files = [f for f in os.listdir(cnst.RTGSFIT_SRC_PATH) if f.endswith('.o')]
    assert not object_files, "Object files were not removed during cleaning."
    # Check constants.c was removed
    constants_c_path = os.path.join(cnst.RTGSFIT_SRC_PATH, 'constants.c')
    assert not os.path.exists(constants_c_path), f"constants.c was not removed: {constants_c_path}"

    # Compile RTGSFIT
    os.chdir(cnst.RTGSFIT_SRC_PATH)
    subprocess.run(
        f"make SHOT={cnst.PULSE_NUM_WRITE} RUN_NAME={cnst.RUN_NAME} DEBUG=1",
        cwd=cnst.RTGSFIT_SRC_PATH,
        check=True,
        shell=True
    )
