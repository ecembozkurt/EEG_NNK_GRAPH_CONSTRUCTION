function X=normalize_ecem(X)
X=(X-min(min(X)))./(max(max(X))-min(min(X)));