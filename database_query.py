import sqlite3
from sqlalchemy import sql

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
        result = cursor.execute(f"SELECT * FROM gene_names WHERE POS = '{search_value}'").fetchall()
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