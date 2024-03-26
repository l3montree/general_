from models import Base

from sqlalchemy import *
from sqlalchemy.orm import *

from pathlib import Path
from typing import *

coinpairs_list_path:Path = "*/database_resources/pairCoins.txt"
coinpair_list_fromfile = None

try:
    with open(coinpairs_list_path,"r") as file:
        coinpairs_list_from_file:List = file.read()
except:
    #creates file if it doesnt exist
    with open(coinpairs_list_path,"x") as file:
        coinpair_list_fromfile = []


class TimeCoinPair(Base):
    __tablename__ = "time_coinpair"

    coinpairs_list = coinpairs_list_from_file

    id = Column(Integer, primary_key=True)
    date_time  = Column(DateTime)
    coinpairs_total = Column(Integer)

    def __setattr__(self,attr):
        self.

    @classmethod
    def append_coinpair_list(cls, coinpair: Union[Integer,List], verbose:bool = False):
        if isinstance(coinpair,List):
            for coinpair_ in coinpair:
                if not coinpair_ in cls.coinpairs_list:
                    cls.coinpairs_list.append(coinpair)
                if verbose:
                    print("Appended {coinpair} to list.\nComplete list = {cls.coinpairs_list}") 
                #update txt file
                with open(coinpairs_list_path, "w") as file:
                    file.write(cls.coinpairs_list)

    @classmethod
    def create_columns(cls,column_list:List):
        new_coinpairs_cols =[]
        for col in column_list:
            if not col in TimeCoinPair.__table__.columns.keys():
                exec(f'{col} = Column(Integer)')
                new_coinpairs_cols.append(col)
        cls.append_coinpair_list(new_coinpairs_cols)

    def __init__(self, datetime:DateTime, coinpairs:Dict[String:Float]):
        date_time = datetime
        coinpairs_total = len(coinpairs)


        for coinpair in coinpairs:
            if coinpair in coinpair_list:









