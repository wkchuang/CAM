import sys
sys.path.insert(0,"/glade/work/wchapman/miniconda3.2/envs/cesmML/lib/python3.11/")
sys.path.insert(1,"/glade/work/wchapman/miniconda3.2/envs/cesmML/lib/python3.11/site-packages")
import numpy as np


def DAMLcnn_run(*args,**kwargs):
    class CNN():
        """
        test this stuff out.
        """

        def __init__(self,args,kwargs):
            super(CNN,self).__init__()
            self.write_quippy_line()
            print("Arguments: ", args)
            print("Keyword arguments: ", kwargs)


        def write_quippy_line(self):
            quippy_line = "What Python Package Env are you Seeing:"
            with open("./describe_env.txt", "w") as file:
                file.write(quippy_line)
                for tt in sys.path:
                    file.write('\n')
                    file.write(tt)
                    file.write('\n')

    CNN(args,kwargs)
    return "Weve done it"
