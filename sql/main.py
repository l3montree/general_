from parse_data import ParseData
from datetime import datetime as _datetime

import argparse
import time

"""
will run parsedata for x seconds
    creates a df from the sql db
        used to create plot of temporal changes in price for each coinpair

"""
datetime_list =[]
curr_datetime: _datetime = None
    
def update_db(timer = 60):
    t_end = time.time()+timer

    while time.time()<t_end:
        curr_data = ParseData()
        curr_data.update_coin_pair()
        curr_datetime = getattr(curr_data, "datetime",None)
        entries = getattr(curr_data,"entries_counter", None)

        if not curr_datetime:
            raise ValueError("no datetime data in current row")

        if not curr_datetime in datetime_list:
            datetime_list.append(curr_datetime)
    print(f'Parse data ran for {timer} s, {entries} appended')


if __name__ == "__main__":
    parser = argparse.add_argument("timer", type = int, help = "Timer (mins) to run the file for", default = None)
    args = parser.parse_args()

    update_db()







    