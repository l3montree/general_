from models import Base

from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy.sql.compiler import IdentifierPreparer

from pathlib import Path
from typing import *
import os.path

import pandas as pd

"""
creates a table with columns = coinpairs, row = datetime, values are prices
    parses coinpair db and updates this table
"""
class TimeCoinPair(Base):
    __tablename__ = "time_coinpair"

    id = Column(Integer, primary_key=True)
    date_time  = Column(DateTime)
    coinpairs_total = Column(Integer)

    coinpairs_list = []

    def __init__(self, datetime:DateTime):
        date_time = datetime
        
        self.coin_pairs_filename = "coin_pairs_list.csv"
        self.coinpairs_list_path:Path = f'database_resources/{self.coin_pairs_filename}'

        if not os.path.isfile(self.coinpairs_list_path):
            open(self.coinpairs_list_path, "w")

        if len(TimeCoinPair.__table__.columns.keys())<4:
            self.load_stored_coinpairs()

        coinpairs_total = len(self.coinpairs_list)

    def load_stored_coinpairs(cls):
        #creates columns for all coinpairs in coin_pairs_list.txt
        stored_coin_pairs = pd.read_csv(cls.coinpairs_list_path)
        stored_coin_pairs = stored_coin_pairs.values.flatten().tolist()

        for coinpair in stored_coin_pairs:
            column = Column(Float, name = coinpair, quote = False)
            TimeCoinPair.__table__.append_column(column)
            cls.coinpairs_list.append(coinpair)


    def get_coinpairs(self):
        return self.coinpairs_list

    

    









