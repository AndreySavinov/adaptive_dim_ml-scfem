function f = right_hand_side(x1, x2)
f = 100*(x1<=0.6).*(x1>=0.4).*(x2<=0.6).*(x2>=0.4);
end