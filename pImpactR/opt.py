import string
import random
import os
from differential_evolution import differential_evolution

def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
    """
    ID = id_generator(size=8)
    input :
      size : (int) character length
    output :
      ID : (string) randomized string
    """
    return ''.join(random.choice(chars) for _ in range(size))
