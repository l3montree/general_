from parse_data import ParseData
from datetime import datetime

datetime_list =[]
curr_datetime: datetime = None

while True:
    curr_data = ParseData()
    curr_data.update_coin_pair()
    curr_datetime = curr_data.get_datetime()

    if not curr_datetime in datetime_list:
        pass
    else:
        pass

    datetime_list.append(curr_data.get_datetime())






    