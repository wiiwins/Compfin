//Solution to C-Exercise 16

//a)
function V_t = BS_Price_DownOut_Call(r, sigma, S_t, T, K, H, t)
    
    a1=(r+sigma^2/2)*(T-t);
    a2=(r-sigma^2/2)*(T-t);
    b=sigma*sqrt(T-t);
    
    if S_t<=H then 
        V_t=0;
    else
        V_t=S_t*(cdfnor("PQ",(log(S_t/K)+a1)/b,0,1)-(H/S_t)^(1+2*r/sigma^2)*cdfnor("PQ",(log(H^2/(K*S_t))+a1)/b,0,1))-exp(-r*(T-t))*K*(cdfnor("PQ",(log(S_t/K)+a2)/b,0,1)-(H/S_t)^(2*r/sigma^2)*cdfnor("PQ",(log(H^2/(K*S_t))+a2)/b,0,1));
    end
    
endfunction
    
//b)
//Setting parameters
r=0.03;
sigma=0.3;
T=1;
H=80;
K=[80,90,100,120];

t=0:0.005:1;
S_t=70:0.5:130;


scf(0)
clf(0)
for i=1:4 // computing prices for strikes = 80,90,100,120
    Price_array=ones(length(S_t),length(t));
    for j=1:length(t)
        for k=1:length(S_t)
            Price_array(k,j)=BS_Price_DownOut_Call(r, sigma, S_t(k), T, K(i), H, t(j));     
        end
    end
       
    subplot(2,2,i)
    surf(t,S_t,Price_array);
    
    //Setting the colormap (thanks to Hyun/Sibanda/Görgen for finding these settings!)
    c=jetcolormap(50)
    c(1,:)=name2rgb('grey')/250 //adding a grey color for V_t around 0
    gcf().color_map = c;
    gce().thickness=0; // setting line width in plot to 0, otherwise plot will be black
    
    label="Strike price: "+string(K(i));
    xtitle(label,"t", "S_t","Value")
    
end
