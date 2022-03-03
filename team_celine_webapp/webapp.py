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
Bootstrap(app)

# index page has form for user's input
@app.route('/', methods=['GET','POST'])
def index():
	if request.method == "POST":
		search_type = request.form.get("search_type") # save search category (snp_id, gene_name, position)
		search_value = request.form.get("search_value") # save search value (text input)
		# go to result page when search button is clicked, passing search type and search value
		return redirect(url_for('search', search_type=search_type, search_value=search_value))
	return render_template('index.html')

# return result data base on search on new url page
@app.route('/search/<search_type>/<search_value>', methods = ['GET'])
def search(search_type, search_value):
	# call function from query script to get requested data
	data = dbq.search_db(search_type, search_value)
	# search type and value needs to be passed since data cannot be passed through routes
	return render_template('result.html', data=data, search_type=search_type, search_value=search_value)

# select statistics and populations form when multiple snps returned (see result.html)
# new url page for statistics results when 'calculate' button is clicked 
@app.route('/stats_pop_select', methods=['POST'])
def stats_pop():
	# get search type and value from search route
	search_type = request.form["search_type"]
	search_value = request.form["search_value"]
	# do another search
	data = dbq.search_db(search_type, search_value) # data from db in list of tuples
	# get selected population and summary statistics as list
	pop_list = request.form.getlist("population")
	stats_list = request.form.getlist("summarystats")
	data_df = dbq.to_df(data) # put data into dataframe for statistics calculation
	# calculate statistics based on user's selection
	stats_df = dbq.calc_stats(data_df, stats_list, pop_list)
	# save statistics dataframe as txt file for the user to download
	stats_df.to_csv('statistics_data.txt', sep=',', index=False, header=True)
	# calculate fst and generate heatmap when multiple populations are selected
	fst_df, heatmap = dbq.calcFst(data_df, pop_list)
	# save fst dataframe as txt file for the user to download
	fst_df.to_csv('fst_data.txt', sep=',', index=True, header=True, index_label="FST")
	# plotting summary statistics
	model_plot = dbq.summary_stats_plot(stats_df, stats_list, pop_list)
	return render_template('stats_pop.html', data=data, search_type=search_type, search_value=search_value,
							pop_list=pop_list, stats_list=stats_list,
							# .to_html() used to automatically show dataframe in a table in html
							# tables=[stats_df.to_html(classes='data', index=False)], fsts=[fst_df.to_html(classes='data', index=True)],
						 	model_plot=model_plot, heatmap=heatmap)

# download txt file
@app.route('/download_statistics')
def download_stats():
	path = "statistics_data.txt" # saved statistics txt file
	return send_file(path, as_attachment=True)

@app.route('/download_fst')
def download_fst():
	path = "fst_data.txt" # saved fst txt file
	return send_file(path, as_attachment=True)

# start the web server
if __name__ == '__main__':
	app.run(debug=True)