library(ggplot2)
ggplot(cars) + 
  aes(x=speed, y=dist) + 
  geom_point() + geom_smooth(method="lm") + labs(title="Speed and Stopping Distances",x="Speed (MPH)", y="Stopping distance (ft)") + theme_bw()

url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
nrow(genes)
colnames(genes)
ncol(genes)
table(genes$State)
table(genes$State)/nrow(genes)*100

p <- ggplot(genes) + 
  aes(x=Condition1, y=Condition2, col=State) + 
  geom_point()

p + scale_color_manual(values=c("blue", "gray", "red"))

library(gapminder)
library(dplyr)
gapminder_2007 <- gapminder %>% filter(year==2007)
ggplot(gapminder_2007) +
  geom_point(aes(x= gdpPercap, y= lifeExp, size=pop), alpha=0.5) + scale_size_area(max_size=10)

ggplot(gapminder) + aes(x=year, y=lifeExp, col=continent) + 
  geom_jitter(width =0.3, alpha=0.4)+ 
  geom_violin(aes(group=year), alpha=0.2, draw_quantiles = 0.5)

gapminder_1957 <- gapminder %>% filter(year==1957)
ggplot(gapminder_1957) + aes(x=gdpPercap, y=lifeExp, col=continent, size=pop) + geom_point(alpha=0.7) +scale_size_area(max_size = 15)

gapminder_1957 <- gapminder %>% filter(year==1957 | year==2007)
ggplot(gapminder_1957) + geom_point(aes(x=gdpPercap, y=lifeExp, color=continent, size = pop), alpha=0.7)+ scale_size_area(max_size = 10) + facet_wrap(~year)

gapminder_top5 <-gapminder %>%
  filter(year==2007) %>%
  arrange(desc(pop)) %>%
  top_n(5, pop)
ggplot(gapminder_top5) + geom_col(aes(x=country, y=pop, fill=lifeExp))
                                       
                                  