from typing import Protocol
import numpy as np

class DataGatheringTool(Protocol):
   def get_data_from_sequence(self, sequence) -> np.ndarray:
       ...
 
   def get_gui_parameters(self) -> dict:
       ...