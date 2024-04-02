from database import Database
from models.CoinPair import CoinPair 
from models.Time_CoinPair import TimeCoinPair

from sqlalchemy import DateTime
from datetime import datetime as _datetime

def temporal_coinpair_changes():
    db = Database()
    with db.db_session() as session:
        time_list = session.query(CoinPair.date_time).distinct()
        sorted(time_list)
     
        for datetime_ in time_list:
            
            datetime_ = datetime_[0]
            time_coinpair = TimeCoinPair(datetime = datetime_)

            coinpairs_each_time = session.query(CoinPair.price,CoinPair.coin_pair).filter(CoinPair.date_time==datetime_).all()
            
            for price,coinpair in coinpairs_each_time:
                setattr(time_coinpair,coinpair,price)
            
            session.add(time_coinpair)

        session.commit()
            

            

        
temporal_coinpair_changes()

"""
methods to do:
 - print out coinpairs vs time for specific period
 - select all coinpairs with highest fluctuations?
    - how to select?
"""


