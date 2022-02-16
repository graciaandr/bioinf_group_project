from cgitb import reset
from re import search
from flask import Flask, render_template, redirect, request, url_for, session
from flask_bootstrap import Bootstrap
import functions_db_query as dbq


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
@app.route('/<search_type>/<search_value>', methods = ['GET', 'POST'])
def search(search_type, search_value):
	data = dbq.search_db(search_type, search_value)
	return render_template('result.html', data=data, search_type=search_type, search_value=search_value)

# select population form
@app.route('/pop_select', methods=['POST'])
def pop_select():
	search_type = request.form["search_type"]
	search_value = request.form["search_value"]
	data = dbq.search_db(search_type, search_value)
	pop_list = request.form.getlist("population")
	# stats_list = request.form.getlist("summarystats")
	filtered_data = dbq.select_pop(pop_list, data)
	# filtered_data = filtered_data.to_numpy()
	return render_template('test.html', data=filtered_data)



# start the web server
if __name__ == '__main__':
	app.run(debug=True, port=5001)