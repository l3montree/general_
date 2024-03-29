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

class ParseData():
    def __init__(self):
        #creates sql database
        self.db = Database()
        self.db.create_database()
        self.price_url = "https://api3.binance.com/api/v3/ticker/price"
        self.data_url = "https://api.binance.com/api/v3/exchangeInfo"
        self.datetime =None
        self.metadata = MetaData()
        self.entries_counter = 0

        self.create_column_iter = 0
    
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

            cp_list =[]

            with self.db.db_session() as session:
                for pair_data in data["symbols"]:
                    _data = _parse_single_coinPair_data(data =pair_data,rename_dict = data_to_rename)
                    _data["price"] = prices_dict.get(_data["coin_pair"]) if True else None
                    
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
                    cp_list.append(cp)
                session.add_all(cp_list)
            self.entries_counter+= len(data["symbols"])
        else:
            print(f"url not parsed")
    
    def print_all_tables(self):
        print(f'Printing all tables and its contents')
        self.metadata = MetaData()
        self.metadata.reflect(bind = self.db.get_engine()) #produces a dictionary of all data in the tables in the engine

        table_names = self.metadata.tables.keys()

        with self.db.db_session(commit = False) as session:
            for table_name in table_names:
                table = table_name
                stmt = select("*").select_from(text(table))
                table_contents = session.execute(stmt)
                print(f'\nTable : {table_name}')
                for row in table_contents:
                    print(row)

    def sql_to_df(self):
        print(f'Converting sql to df')
        self.metadata = MetaData()
        self.metadata.reflect(bind = self.db.get_engine()) #produces a dictionary of all data in the tables in the engine

        table_names = self.metadata.tables.keys()
        
        file_path:Path = "*/database_resources/df_files/"
        df_files_in_dir = glob.glob("*/database_resources/df_files/*")
        file_ext =".xlsx"

        """
        check if db is present
            if not create
            else append df

        """

        with self.db.db_session(commit = False) as session:
            for table_name in table_names:
                table = table_name
                stmt = select("*").select_from(text(table))
                df = pd.read_sql(stmt, session.bind)

                if table_name in df_files_in_dir:
                    df2 = pd.read(f'{file_path}{table_name}{file_ext}')
                    df3 = pd.concat([df2,df])




                df.to_excel(f'{table_name}{file_ext}')



if __name__ == "__main__":
    parse_data = ParseData()





