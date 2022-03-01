from flask import Flask, render_template, redirect, url_for, send_file, request
from flask_bootstrap import Bootstrap
import pandas as pd
import functions_db_query as dbq
import matplotlib.pyplot as plt 

# a non-interactive backend that can only write to files
# used on Linux, if Matplotlib cannot connect to either an X display or a Wayland display
plt.switch_backend('Agg')

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
	return render_template('result.html', data=data, search_type=search_type, search_value=search_value)

# select statistics and populations form and return calculate statistics
@app.route('/stats_pop_select', methods=['POST'])
def stats_pop():
	search_type = request.form["search_type"]
	search_value = request.form["search_value"]
	data = dbq.search_db(search_type, search_value) # data from db in list of tuples
	pop_list = request.form.getlist("population")
	stats_list = request.form.getlist("summarystats")
	data_df = dbq.to_df(data) # put data into dataframe for calculation
	# calculate statistics based on user's selection
	stats_df = dbq.calc_stats(data_df, stats_list, pop_list)
	stats_df.to_csv('statistics_data.txt', sep=',', index=False, header=True)
	# calculate fst and generate heatmap when multiple populations are selected
	fst_df, heatmap = dbq.calcFst(data_df, pop_list)
	fst_df.to_csv('fst_data.txt', sep=',', index=True, header=True, index_label="FST")
	# plotting summary statistics
	model_plot = dbq.summary_stats_plot(stats_df, stats_list, pop_list)
	return render_template('stats_pop.html', data=data, search_type=search_type, search_value=search_value,
							pop_list=pop_list, stats_list=stats_list,
							tables=[stats_df.to_html(classes='data', index=False)], fsts=[fst_df.to_html(classes='data', index=True)],
						 	model_plot=model_plot, heatmap=heatmap)

# download txt file
@app.route('/download_statistics')
def download_stats():
	path = "statistics_data.txt"
	return send_file(path, as_attachment=True)
@app.route('/download_fst')
def download_fst():
	path = "fst_data.txt"
	return send_file(path, as_attachment=True)

# start the web server
if __name__ == '__main__':
	app.run(debug=True, port=5001)