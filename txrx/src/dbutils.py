import pymongo
import pandas as pds
import numpy as np;
import qgrid
from bson.objectid import ObjectId

_myclient = pymongo.MongoClient("mongodb://localhost:27017/")
db_exp = _myclient["test"]['experiments']
db_src = _myclient['test']['source_code']
as_table = lambda x : pds.DataFrame(list(x));


from hashlib import md5
from datetime import datetime
from copy import deepcopy
import inspect

def get_wf(x):
	data = np.fromfile(f'data/{x}.bin', dtype='float32')
	dbe = list(db_exp.find({'_id' : ObjectId(x)}));
	if len(dbe) and 'rx_ch' in dbe[0]:
		data = data.reshape([-1, len(dbe[0]["rx_ch"])]);
	return data
	
def try_save_src(fn):
	fn_src = inspect.getsource(fn);
	fn_hash = md5(fn_src.encode('utf-8')).hexdigest();
	if db_src.count_documents({'fn_hash' : fn_hash}) == 0: 
		rec = db_src.insert_one({'fn_hash' : fn_hash, 'datetime' : str(datetime.now()), 'fn_name' : fn.__qualname__, 'fn_src' : fn_src })			
		print(f'saved new or revised function:  {fn.__qualname__} inserted as {rec.inserted_id} ');
def save_src_decor(fn):
	def wrapper(*args, **kwargs):
		try_save_src(fn);
		return fn(*args, **kwargs);
	return wrapper;
def np_change_dtype(di):
	di2 = {}
	for k,v in di.items():
		if getattr(v, "item", None):
			di2[k] = v.item();
		else:
			di2[k] = v;
	return di2;

def save_result_decor(fn):
	def wrapper(**kwargs):# no *args, user must not use positional argument (use keyword argument instead)
		info = deepcopy(kwargs);
		try_save_src(fn);
		info['call_fn'] = fn.__qualname__;
		insert_res = db_exp.insert_one(info);
		fid = insert_res.inserted_id;	   
		data = fn(**kwargs);
		now = str(datetime.now()); # note the datetime here should >= datetime of saved source, so we can $max the source database
		# to find the function related to certain record.
		db_exp.update_one({'_id' : ObjectId(fid)}, {'$set': {'datetime' : now}});
		if (type(data) is np.ndarray):
			data.tofile(f'data/{fid}.bin');
		elif type(data) is dict: # for processing results (like BER, snr etc, insert the result in the record)
			db_exp.update_one({'_id' : ObjectId(fid)}, {'$set': np_change_dtype(data)});
		else:
			raise ValueError("Unexpected return value, require np.ndarray(save to file) or dict(save to database)");
		print(f'experiment data saved _id={fid}, date={now}')
		return data;
	return wrapper;


#example of experiment function:
#any argument will be saved, in order to document any relevant variables of experiment
#return a numpy array, will be saved as (hex id).bin, whose id is the _id field of recorded metainfo
#because the wr
@save_result_decor
def test_db(arg1, rx_distance, rotation_angle, tx_distance, temperature, wind_speed):
	print(f'Hello, World! arg1 is {arg1}');
	return np.zeros((2,2));
