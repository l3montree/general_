"""
Tutorial:
links:
- https://docs.sqlalchemy.org/en/20/tutorial/engine.html

google doc: https://docs.google.com/document/d/17KuoHsSyw7nVf0xZZhqaFAz6roCo6NnBN0yYpd7UZOQ/edit 


"""
from sqlalchemy import create_engine, text
engine = create_engine("sqlite+pysqlite:///:memory:", echo=True) #creates in memory (local) database

with engine.connect() as conn:
    conn.execute(text("CREATE TABLE some_table (x_ int, y_ int)")) #converts strings into sql commands
    # CREATE {database_type} {database_name} ({columns_name} {column_dtype})
    # {ACTION} {params}
    conn.execute(text("INSERT INTO some_table (x_, y_) VALUES (:x_, :y)"), [{"x_": 1, "y": 1}, {"x_": 2, "y": 4}] ,)
    # {ACTIOn: INSERT} {directional: INTO} {database_name} {columns} {VALUES --> input params} {param_names} \n {param_vals corresp to param_names}
    conn.commit()
    results = conn.execute(text("SELECT x_,y_ FROM some_table"))
    for i in results:
        print(i)

#QUERY