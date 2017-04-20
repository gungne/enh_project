find . -size +10M | sed 's|^\./||g' | cat >> .gitignore
find . -type f -name "*.fa"| sed 's|^\./||g' | cat >> .gitignore
git add .
git commit -m 'update'
git push