
des1 <- data.frame(A = factor(c(1,2,3)))
print(suppressWarnings(generate_main_effect_contrast(des1, "A", cellcount=rep(1,nrow(des1)))))

des2 <- data.frame(A = factor(c(1,2,3,4)))
print(suppressWarnings(generate_main_effect_contrast(des2, "A", cellcount=rep(1,nrow(des2)))))

des3 <- expand.grid(A =  c("1","2","3","4"), B=c("a", "b", "c"))
print(des3)
print(suppressWarnings(generate_main_effect_contrast(des3, "A", cellcount=rep(1,nrow(des3)))))
print(suppressWarnings(generate_main_effect_contrast(des3, "B", cellcount=rep(1,nrow(des3)))))

des3 <- expand.grid(A =  c("1","2","3","4"), B=c("a", "b", "c"))
print(des3)
print(generate_interaction_contrast(des3, c("A", "B")))