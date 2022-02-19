from cgitb import reset
from os import stat
from re import search
from flask import Flask, render_template, redirect, request, url_for, session
from flask_bootstrap import Bootstrap
import pandas as pd
import functions_db_query as dbq


# create a flask application object
app = Flask(__name__)

# flask_wft requires encryption key
app.config['SECRET_KEY'] = 'chromosome22'
Bootstrap(app)


# search form in home page
@app.route('/', methods=['GET','POST'])
def index():
	if request.method == "POST":
		search_type = request.form.get("search_type")
		search_value = request.form.get("search_value")
		return redirect(url_for('search', search_type=search_type, search_value=search_value))
	return render_template('index.html')

# return result data base on search
@app.route('/search/<search_type>/<search_value>', methods = ['GET'])
def search(search_type, search_value):
	data = dbq.search_db(search_type, search_value)
	# data_df = dbq.to_df(data)
	# df_cols = data_df.columns
	return render_template('result.html', data=data, search_type=search_type, search_value=search_value)

# select statistics and populations form
@app.route('/stats_pop_select', methods=['POST'])
def stats_pop():
	search_type = request.form["search_type"]
	search_value = request.form["search_value"]
	data = dbq.search_db(search_type, search_value) # data from db in list of tuples
	pop_list = request.form.getlist("population")
	stats_list = request.form.getlist("summarystats")
	data_df = dbq.to_df(data) # put data into dataframe for stats

	if 'shannon' in stats_list:
		shannon_df = pd.DataFrame({'pop' : pop_list}) # make df for shannon value per population
		shannon_list = []
		for pop in pop_list: # loop to get shannon per pop
			shannon = dbq.calcShannonDiv(data_df, pop)
			shannon_list.append(shannon)
		shannon_df["shannon"] = shannon_list # input shannon values into df
		shannon = shannon_df.to_numpy() # df to list of tuples for html
	
	elif 'tajima' in stats_list:
		for pop in pop_list:
			tajima = dbq.calcTajimaD(data_df, pop)
	
	# elif 'hetero' in stats_list:

	else:
		print ("Choose statistics to calculate")



	return render_template('stats_pop.html', data=data, shannon=shannon)


# start the web server
if __name__ == '__main__':
	app.run(debug=True, port=5001)