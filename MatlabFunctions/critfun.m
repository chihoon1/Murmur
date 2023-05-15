function dev = critfun(X,Y)
model = fitglm(X,Y,'Distribution','binomial');
dev = model.Deviance;
end     