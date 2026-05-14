from enum import Enum


class Origin(Enum):
    THAI = 0
    MIGRANT = 1

class Route(Enum):
    UNKNOWN = 0
    VERTICAL = 1
    HORIZONTAL_HH = 2
    HORIZONTAL_COM = 3