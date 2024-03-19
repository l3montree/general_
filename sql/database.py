import sqlalchemy.engine.url as url
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker,scoped_session
from models import *

from contextlib import contextmanager

class Database:
    def __init__(self):
        self.uri  = "sqlite:///database_resources/tutorial.db" #local db name
        self.engine = create_engine(self.uri)
        self.SessionMaker = sessionmaker(bind = self.engine)
        
    def get_metadata(self):
        return Base.metadata
    
    def get_engine(self):
        return self.engine
        
    def create_database(self):
        Base.metadata.create_all(self.engine) #creates database
    
    @contextmanager
    def db_session(self, commit = True):
        session:Session = scoped_session(self.SessionMaker)
        yield session #code exits prematurely to complete all code after the function call, then returns here to complete the remaining code in this function
        
        if commit:
            session.commit()
        session.close()


        

        