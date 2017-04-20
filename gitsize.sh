 > .gitignore
find . -size +10M | sed 's|^\./||g' | cat >> .gitignore
# sort .gitignore|uniq >.gitignore
git add .
git commit -m 'update'
git push