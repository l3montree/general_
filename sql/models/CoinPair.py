from models.base import Base

from sqlalchemy import Column, ForeignKey, Integer, String, Float, CheckConstraint, DateTime, Boolean, event
from sqlalchemy.orm import relationship, Session

"""
SQL table class stores: coins,coinpairs, price etc
    stores raw data that will be analysed later
"""

class CoinPair(Base):
    __tablename__ ="coin_pair"

    id = Column(Integer,primary_key=True)
    coin_pair = Column(String)
    from_coin = Column(String)
    to_coin = Column(String)
    price =Column(Float)
    status = Column(String)
    enabled = Column(Boolean)
    
    date_time = Column(DateTime)

    __table_args__ =(
        CheckConstraint('from_coin <> to_coin', name ="different_coins_constraint"),)#ensures to_coin != from_coin

    def __init__(self, coin_pair,from_coin,to_coin,date_time,status,enabled,price =None):
        self.coin_pair = coin_pair
        self.to_coin =to_coin
        self.from_coin  = from_coin
        self.price = price
        self.status = status
        self.date_time = date_time
        self.enabled = enabled

    
