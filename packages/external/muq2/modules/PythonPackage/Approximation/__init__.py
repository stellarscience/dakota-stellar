import sys
from ..pymuqModeling import ModPiece
from ..pymuqApproximation import *

sys.modules[__name__] = sys.modules['muq.pymuqApproximation']
