import sqlite3
from sqlite3 import Error


def create_db(db_file):
    """ create a database connection to a SQLite database """
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except Error as e:
        print(e)
    finally:
        conn.close()


def connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print "No existing database, please create new database"

    return None

def create_table(conn,tablename,columns):
    template = "CREATE TABLE IF NOT EXISTS _table_name_ ( _column_desc_ );"

    sql = template.replace("_table_name_",tablename)
    sql = sql.replace("_column_desc_", ", ".join(columns))

    try:
        c = conn.cursor()
        c.execute(sql)
    except Error as e:
        print(e)



if __name__ == '__main__':
    geneIDDatabase = readData()
    print data["ENSG00000005073"]