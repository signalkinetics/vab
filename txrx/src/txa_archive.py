def constel(data, text = None,wid=500):
    if (text == None) :
        text = list(map(str, np.arange(data.shape[0]))); 
    fig=go.Figure();
    fig.add_trace(go.Scatter(
        x=np.real(data),
        y=np.imag(data),
        text=text,
        mode="markers+text",
        textposition="bottom center"
    ))
    fig.update_layout(autosize=False, width=wid,    height=wid,)
    fig.update_yaxes(scaleanchor='x')
    return fig;


def txa_goodness(offset):
    pair = dat1[offset : offset + spb * 6];
    
    window = chebwin(spb, 150).T;
    W = lambda x : x @ window;
    a = W(pair[int(spb*.5):int(spb*1.5)]);
    ab = W(pair[int(spb*2.5):int(spb*3.5)]);
    b = W(pair[int(spb*4.5):int(spb*5.5)]);
    normer = np.sqrt(np.abs(a) * np.abs(b));
    return {'a' : a, 'b' : b, 'ab' : ab, 'offset': offset,
        'amplitude' : -np.abs(np.abs(a) - np.abs(b)) / normer,
            'phase' : np.real(a * np.conj(b)) / (np.abs(a) * np.abs(b))-1, 
            'linear' : -np.abs(ab - (a+b))/normer}

df = pds.DataFrame([txa_goodness(spb * i + 1700 * j +4) for i in [4,8,12,16,20,26] for j in range(1,20)])
import plotly.express as px
fig = px.ecdf(df, x=['amplitude','phase', 'linear'])
fig.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),width=400,
);fig.update_xaxes(range=[-2,0])
fig.show(renderer='notebook')
data_dir = 'data/20220729/'
for rot in ['0', 'm30','p30']:
    segment=np.fromfile(f'{data_dir}_50V_3m_rot_{rot}.bin', dtype='float32');
    sm = mixer(Fc, Fs, segment, 400).reshape((-1,)); pc(sm);
    val = []; label = []; si = 0;
    for r in [1,2,3,4]:
        for k in [1]:
            for i in [1,2]:
                si = i - 1 + k * 3 + r * 12;
                #d = fft(window * segment[si * 50000 + 12500:si * 50000 + 37500]); val.append(d[2500]); 
                val.append(sm[50+si*100])
                label.append(f'r{r}_k{k}_i{i}')
                
    constel(val,label,1000).show(render='notebook')
    
data_dir = 'data/20220716/'
for dist in range(1,6):
    segment = np.fromfile(f'{data_dir}_25V_{dist}m_rot0.bin', dtype='float32')
    window=chebwin(25000, 100);
    sm = mixer(Fc, Fs, segment, 1000).reshape((-1,)); pc(sm);
    for j in [[15]+list(range(1,15)), range(15,30)]:
        val = [0]*16
        if(j[-1]==29):
            continue;
        for i in range(1,16):
            si = j[i - 1];
            d = fft(window * segment[si * 50000 + 12500:si * 50000 + 37500])
            val[i]= d[2500];
        constel(np.array(val)).show(render='notebook')
