function [l,u] = piecewise_linear_bounds_tanh(W,b,X,epsilon)
num_layers = length(W)-1;
[ll,uu]=first_pre_activate(W{1},b{1},X,epsilon);
l=[];
u=[];
l=[l;ll];
u=[u;uu];
leng=length(ll);
for i=1:num_layers-1
    [ll,uu]=next_layer_prebound(W{i+1},b{i+1},ll,uu);
    l=[l;ll];
    u=[u;uu];
end


    function  [l_first, u_first] = first_pre_activate(W,b,X,epsilon)
        l_first=-abs(W)*epsilon  +  W*X  +  b;
        u_first= abs(W)*epsilon  +  W*X  +  b;
    end



    function   y_prim = tanh_diff(x)
        epsi=0.0001;
        y_prim= (tanh(x+epsi)-tanh(x-epsi)) /  (2*epsi)  ;
    end




    function   [l_next, u_next] = next_layer_prebound(W,b,l_pre,u_pre)
        leng=length(l_pre);
        alph_upp=[];
        alph_low=[];
        alphabet_upp=[];
        alphabet_low=[];
        for i=1:leng
            if l_pre(i,1)>0
                d=0.5*(  l_pre(i,1) +   u_pre(i,1)   );
                alph_upp=[alph_upp;   tanh_diff(d)];
                alphabet_upp=[alphabet_upp;  tanh(d)-tanh_diff(d)*d];
                ss=(  tanh(u_pre(i,1))-tanh(l_pre(i,1))  ) / (    u_pre(i,1)-l_pre(i,1)   );
                alph_low=[ alph_low ; ss ];
                alphabet_low=[ alphabet_low ; tanh(l_pre(i,1))-ss*l_pre(i,1)  ];
            elseif u_pre(i,1)<0
                ss=(tanh(u_pre(i,1))-tanh(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
                alph_upp=[alph_upp ; ss ];
                alphabet_upp=[alphabet_upp ; tanh(l_pre(i,1))-ss*l_pre(i,1)   ];
                d=0.5*(l_pre(i,1)+u_pre(i,1));
                alph_low=[alph_low ; tanh_diff(d) ];
                alphabet_low=[alphabet_low ; tanh(d)-tanh_diff(d)*d ];
            else
                x0=0.01;
                fun=@(x)d_find(x,l_pre(i,1));
                option=optimset('MaxFunEvals',500000);
                option=optimset(option,'MaxIter',500);
                option=optimset(option,'disp','iter','LargeScale','off','TolFun',0.0001);
                [d,fval]=fsolve(fun,x0,option);
                clc
                if (d<0 || d>u_pre(i,1))
                    ss=(tanh(u_pre(i,1))-tanh(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
                    alph_upp=[alph_upp  ;  ss   ];
                    alphabet_upp=[alphabet_upp  ; tanh(l_pre(i,1))-ss*l_pre(i,1)  ];
                else
                    alph_upp=[alph_upp  ; tanh_diff(d)];
                    alphabet_upp=[alphabet_upp  ; tanh(l_pre(i,1))-tanh_diff(d)*l_pre(i,1) ];
                end
                
                x0=0.01;
                fun=@(x)d_find(x,u_pre(i,1));
                option=optimset('MaxFunEvals',500000);
                option=optimset(option,'MaxIter',500);
                option=optimset(option,'disp','iter','LargeScale','off','TolFun',0.0001);
                [d,fval]=fsolve(fun,x0,option);
                clc
                if (d>0 || d<l_pre(i,1))
                    ss=(tanh(u_pre(i,1))-tanh(l_pre(i,1)))/(u_pre(i,1)-l_pre(i,1));
                    alph_low=[alph_low ; ss  ];
                    alphabet_low=[alphabet_low ; tanh(u_pre(i,1))-ss*u_pre(i,1) ];
                else
                    alph_low=[alph_low ; tanh_diff(d)];
                    alphabet_low=[alphabet_low ; tanh(u_pre(i,1))-tanh_diff(d)*u_pre(i,1)   ];
                end
            end
        end
        
        
        
        
        W_pos=0.5*(W+abs(W));
        W_neg=0.5*(W-abs(W));
        mid=0.5*(l_pre+u_pre);
        dif=0.5*(u_pre-l_pre);
        cc1=  W_pos.*alph_low.'   +    W_neg.*alph_upp.';
        dd1=  W_pos*alphabet_low  +    W_neg*alphabet_upp + b;
        l_next=-abs(cc1)*dif   + cc1*mid + dd1 ;
    
    
    
        cc2=  W_neg.*alph_low.'    +    W_pos.*alph_upp.';
        dd2=  W_neg*alphabet_low   +    W_pos*alphabet_upp + b;
        u_next= abs(cc2)*dif   +  cc2*mid + dd2;
        

    end



    function Eq = d_find(x,data)
        
        Eq=((tanh(x)-tanh(data))/(x-data))-tanh_diff(x);
        
    end

end
