from flask import Flask, render_template, redirect, request, url_for
from flask_bootstrap import Bootstrap
import database_query as dbq


app = Flask(__name__)

# flask_wft requires encryption key
app.config['SECRET_KEY'] = 'anythingreally'
Bootstrap(app)

# put form in route function
@app.route('/', methods=['GET','POST'])
def index():
	snp_id = None
	if request.method == "POST":
		search_type = request.form.get("search_type")
		search_value = request.form.get("search_value")
		data = dbq.search_db(search_type, search_value)
		return render_template('index.html', data=data)
	return render_template('index.html')


# start the web server
if __name__ == '__main__':
	app.run(debug=True, port=5001)