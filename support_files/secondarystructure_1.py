# perform secondary structure
# transfer to GUI
# annotation for components on secondary structure with %similarity against ref

## project_alignate/external_tools/psipred/BLAST+
## project_alignate/external_tools/psipred/blastdb

from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout
from PySide6.QtCore import Signal


class build_secondary_structure(QWidget):

    #Signal to send data back to alignate.py
    result_ready = Signal(str)

    def __init__(self):
        super().__init__()

        # 1 Get sequences to run
        dir_path = os.path.dirname(os.path.abspath(__file__))           # full/absolute path to current dir (full path includes filename)
        eg_onesamplefasta = os.path.join(dir_path, 'example_datasets', 'onesample.fasta')
        eg_twosamplefasta = os.path.join(dir_path, 'example_datasets', 'twosamples.fasta')

        # 2 Run psipred
        psipred_path = os.path.join(dir_path, 'external_tools', 'psipred', 'runpsipred_single')
        try:
            subprocess.run([psipred_path, eg_onesamplefasta], check=True)
            print('✅ PSIPRED run successful')
        except subprocess.CalledProcessError as e:
            print(f'❌ PSIPRED failed: {e}')

        # 3 Display on GUI
        print('secondarystructure done')



        self.result_ready.emit(self.result)

        