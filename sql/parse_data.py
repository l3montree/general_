import requests
from database import Database

from typing import Dict,Union
from sqlalchemy.orm import Session
from sqlalchemy import event, MetaData, select, text, Column, types, DateTime

from models import *

from datetime import datetime as _datetime

from time import process_time


import pandas as pd
from pathlib import Path
import glob
import os

import xlsxwriter

class ParseData():
    def __init__(self,verbose = False):
        #creates sql database
        self.db = Database()
        self.db.create_database()
        self.price_url = "https://api3.binance.com/api/v3/ticker/price"
        self.data_url = "https://api.binance.com/api/v3/exchangeInfo"
        self.datetime =None
        self.metadata = MetaData()
        self.entries_counter = 0

        self.create_column_iter = 0

        self.verbose = verbose

        self.coinpairs_list =[]
        self.coinpairs_list_path_txt = "database_resources/coin_pairs_list.txt"
        self.coinpairs_list_path_csv = "database_resources/coin_pairs_list.csv"
    
    def get_datetime(self):
        return self.datetime
        
    def update_coin_pair(self):
        def _parse_single_coinPair_data(data:Dict[str,str],rename_dict: Dict[str,str]):
            parsed_data: Dict[str:str] = {}
            for key,value in rename_dict.items():
                if key in data:
                    parsed_data[value] = data[key]
            
            return parsed_data
        

        # Make the HTTP GET request
        data_response = requests.get(self.data_url)
        price_response = requests.get(self.price_url)

        # updates coinpairs:
        # Check if the request was successful (status code 200)
        if data_response.status_code == 200 & price_response.status_code==200:
            # Parse the JSON response
            data = data_response.json()
            prices =price_response.json()
            prices_dict ={}


            # Dict[symbol:prices]
            for i in prices:
                key = i["symbol"]
                value = i["price"]
                if not key in prices_dict.keys():
                    prices_dict[key] = value
                else:
                    raise ValueError(f"duplicate coin_pairs:{key}")

            # Extract the price from the response
            data_to_rename = {"symbol":"coin_pair","baseAsset":"from_coin","quoteAsset":"to_coin","status":"status"}
            
            self.datetime =_datetime.now()

            cp_list =[] #stores coinpair: sql instances
            coinpairs = [] #stores all coinpairs:str passed to cp_list

            with self.db.db_session() as session:
                for pair_data in data["symbols"]:
                    _data = _parse_single_coinPair_data(data =pair_data,rename_dict = data_to_rename)
                    _data["price"] = prices_dict.get(_data["coin_pair"]) if True else None
                    
                    if not _data["coin_pair"] in self.coinpairs_list:
                        self.coinpairs_list.append(_data["coin_pair"])
                    #enabled column
                    if _data["status"].lower() == "trading":
                        _enabled = True
                    else:
                        _enabled = False
                    
                    #uploads all coinpairs
                    
                    cp =CoinPair(
                        coin_pair = _data["coin_pair"], to_coin = _data["to_coin"], from_coin = _data["from_coin"],\
                        status = _data["status"], date_time = self.datetime, enabled = _enabled, price = _data["price"]
                    )
                    cp_list.append(cp) #all coin pairs appended here
                    
                    coinpairs.append(_data["coin_pair"])
                    self.entries_counter+=1
                
                session.add_all(cp_list)

            #update coinpair_list.txt
            coinpairs = list(set(coinpairs))
            self.update_coinpairs_csv(coinpairs)
        else:
            print(f"url not parsed")


    def update_time_coinpair(self):

        with self.db.db_session() as session:
            time_list = session.query(CoinPair.date_time).distinct()
            sorted(time_list)

            out = 0

            tc = time_coinpair = TimeCoinPair(datetime = _datetime.now())
            session.add(tc)
            session.commit()

        
            for datetime_ in time_list:
                out +=1
                
                datetime_ = datetime_[0]
                time_coinpair = TimeCoinPair(datetime = datetime_)

                coinpairs_each_time = session.query(CoinPair.price,CoinPair.coin_pair).filter(CoinPair.date_time==datetime_).all()
                
                for price,coinpair in coinpairs_each_time:
                    setattr(time_coinpair,coinpair,price)
                session.add(time_coinpair)
                session.commit()

            


    #misc methods

    def update_coinpairs_csv(self, coinpairs_list:list):

        if not os.path.isfile(self.coinpairs_list_path_csv):
            open(self.coinpairs_list_path_csv,"w")
            
        
        coinpairs_list = pd.Series(coinpairs_list)
        
        # retrieving stored csv 
        try:
            stored_coinpairs = pd.read_csv(self.update_coinpairs_csv) 
        except:
            print(f'parsedata.update_coinpairs_csv: pd.read_csv(path) has an empty csv')
            stored_coinpairs = pd.Series([]) #creates an empty pd.series 


        if not len(stored_coinpairs):
            print("coinpair_list.csv empty, initiatlising now")
        else:
            coinpairs_list.append(stored_coinpairs)

        if len(coinpairs_list):
            coinpairs_list.sort_values(inplace = True)
            coinpairs_list.to_csv(self.coinpairs_list_path_csv, index = False)

        else:
            raise ValueError("Parse_data.update_coin_pairs_txt received a null list + no stored coin pairs")

    def sql_to_excel(self):
        print(f'Converting sql to df')
        self.metadata = MetaData()
        self.metadata.reflect(bind = self.db.get_engine()) #produces a dictionary of all data in the tables in the engine

        table_names = self.metadata.tables.keys()
        
        file_path = Path("database_resources/df_files")
        df_files_in_dir = glob.glob(os.path.join(file_path,"*"))
        file_ext =".xlsx"

        with self.db.db_session(commit = False) as session:
            for table_name in table_names:
                table = table_name
                stmt = select("*").select_from(text(table))
                df = pd.read_sql(stmt, session.bind)

                df_filepath = os.path.join(file_path,f'{table_name}{file_ext}')
                df.to_excel(df_filepath, index = False)
                


if __name__ == "__main__":
    parse_data = ParseData()





