install.packages("installr")
installr::updateR()

# 创建一个示例数据框
df <- data.frame(
    name = c("Alice", "Bob", "Charlie"),
    age = c(25, 30, 35),
    score = c(85, 90, 95)
)

# 使用 View 函数查看数据框
View(df)
