from database import Database

def check_fluctuations(timeList):
    db = Database()
    with db.db_session() as session:
        pass


"""
methods to do:
 - print out coinpairs vs time for specific period
 - select all coinpairs with highest fluctuations?
    - how to select?
"""

