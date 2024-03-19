from models import *

from datetime import datetime
from database import Database

db = Database()
db.create_database()

time_now = datetime.now()
with db.db_session() as session:
    t = Time(date_time=time_now)
    

