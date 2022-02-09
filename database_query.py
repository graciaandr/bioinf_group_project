import sqlite3, csv
from sqlalchemy import sql

# insert values into table in db
def insert_table(csv_file, table_name):
    with sqlite3.connect("chr22.db") as con:
        cursor = con.cursor()
        with open (csv_file, 'r') as i:
            reader = csv.reader(i)
            columns = next(reader) 
            query = "INSERT INTO " + table_name + " ({0}) VALUES ({1})"
            query = query.format(','.join(columns), ','.join('?' * len(columns)))
            for data in reader:
                cursor.execute(query, data)
            con.commit()



############### query functions ################

def use_id(search_value):
    with sqlite3.connect("chr22.db") as connection:
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM gene_names WHERE ID = '{search_value}'").fetchall()
    return result

def use_gene(search_value):
    with sqlite3.connect("chr22.db") as connection:
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM gene_names WHERE GENE = '{search_value}'").fetchall()
    return result

def use_pos(search_value):
    with sqlite3.connect("chr22.db") as connection:
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM GBR_genotype WHERE POS = '{search_value}'").fetchall()
    return result

def search_db(search_type, search_value):
    if search_type == "snp_id":
        return use_id(search_value)
    elif search_type == "gene_name":
        return use_gene(search_value)
    elif search_type == "genomic_coordinate":
        return use_pos(search_value)
    else:
        return "404: Not Found"