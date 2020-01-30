library(tidyverse)
library(microbenchmark)

source("utils.R")
source("utils_old.R")


###############################################################
###############################################################
# Test speedup for the computation of the Scheffé model matrix
###############################################################
###############################################################

n_runs = sample(5:200, 1)
q = sample(3:10, 1)
X = create_random_initial_design(n_runs, q)
# X = create_random_initial_design(5, 3)
r_ord2 = get_scheffe_order_2(X)
c_ord2 = getScheffeOrder2(X)

print(n_runs)
print(q)
sum(r_ord2 - c_ord2)






n_runs = sample(5:200, 1)
q = sample(3:10, 1)
X = create_random_initial_design(n_runs, q)
# X = create_random_initial_design(5, 3)
r_ord3 = get_scheffe_order_3(X)
c_ord3 = getScheffeOrder3(X)

print(n_runs)
print(q)
sum(r_ord3 - c_ord3)





r_ord3 = get_scheffe(X, 3)
c_ord3 = getScheffe(X, 3)
sum(c_ord3 - r_ord3)




n_runs = sample(5:200, 1)
q = sample(3:10, 1)
test_order_3 = function(n_runs = 10, q = 3){
  
  X = create_random_initial_design(n_runs, q)
  r_ord3 = get_scheffe_order_3(X)
  c_ord3 = getScheffeOrder3(X)
}



q = 3
results_q_3 = lapply(seq(100, 5000, by = 100), function(n_runs){
  print(n_runs)
  time_res = microbenchmark(
    {
      X = create_random_initial_design(n_runs, q, seed = n_runs)
      
      r_ord3 = get_scheffe_order_3(X)
    },
    {
      X = create_random_initial_design(n_runs, q, seed = n_runs)
      
      c_ord3 = getScheffeOrder3(X)
    },
    times = 50
  )
  time_res %>% 
    as.data.frame() %>% 
    group_by(expr) %>% 
    summarize(
      mean = mean(time),
      median = median(time),
      p25 = as.numeric(quantile(time, 0.25)),
      p75 = as.numeric(quantile(time, 0.75))) %>% 
    mutate(n_runs = n_runs)
})


results_q_3 %>% 
  bind_rows() %>% 
  ggplot(aes(x = n_runs, y = median, 
             ymin = p25,
             ymax = p75,
             group = expr,
             fill = expr)) +
  geom_line(aes(color = expr)) +
  geom_ribbon(alpha = 0.3) +
  theme_bw()



results_q_3 %>% 
  bind_rows() %>% 
  select(median, expr, n_runs) %>% 
  pivot_wider(values_from = c(median), 
              names_from = expr) %>% 
  mutate(speedup = r_ord3/c_ord3) %>% 
  ggplot(aes(x = n_runs, y = speedup)) +
  geom_line() +
  theme_bw()






# results_all = lapply(seq(10, 5000, by = 100), function(n_runs){
results_all = lapply(seq(10, 400, by = 100), function(n_runs){
  lapply(seq(5, 20, by = 5), function(q){
    cat("q =", q, ", n_runs =", n_runs, "\n")
    time_res = microbenchmark(
      R_expr = {
        X = create_random_initial_design(n_runs, q, seed = n_runs)
        r_ord3 = get_scheffe_order_3(X)
      },
      cpp_xpr = {
        X = create_random_initial_design(n_runs, q, seed = n_runs)
        c_ord3 = getScheffeOrder3(X)
      },
      times = 10
    )
    time_res %>% 
      as.data.frame() %>% 
      group_by(expr) %>% 
      summarize(
        mean = mean(time),
        median = median(time),
        p25 = as.numeric(quantile(time, 0.25)),
        p75 = as.numeric(quantile(time, 0.75))) %>% 
      mutate(n_runs = n_runs,
             q = q)
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()




df_plot = results_all %>% 
  select(median, expr, n_runs, q) %>% 
  pivot_wider(values_from = c(median), 
              names_from = expr) %>% 
  mutate(median_speedup = R_expr/cpp_xpr,
         q = as.factor(q)) %>% 
  select(n_runs, q, median_speedup) %>% 
  inner_join(
    results_all %>% 
      select(mean, expr, n_runs, q) %>% 
      pivot_wider(values_from = c(mean), 
                  names_from = expr) %>% 
      mutate(mean_speedup = R_expr/cpp_xpr,
             q = as.factor(q)) %>% 
      select(n_runs, q, mean_speedup)
  )


df_plot %>% 
  ggplot(aes(x = n_runs, y = median_speedup, color = q)) +
  geom_line() +
  theme_bw()


df_plot %>% 
  ggplot(aes(x = n_runs, y = median_speedup, color = q)) +
  geom_line() +
  facet_wrap(~q, scales = "free_y") +
  theme_bw()



df_plot %>% 
  ggplot(aes(x = n_runs, y = mean_speedup, color = q)) +
  geom_line() +
  theme_bw()


df_plot %>% 
  ggplot(aes(x = n_runs, y = mean_speedup, color = q)) +
  geom_line() +
  facet_wrap(~q, scales = "free_y") +
  theme_bw() 



###############################################################
###############################################################
# Test speedup for the whole algorithm using finite approximations
###############################################################
###############################################################

# Takes like half an hour.
# 1.2x speedup
time_res_1 = microbenchmark(
  e_C = {
    res_alg_order_C = coord_ex_mixt(
      20, q = 15, 
      order = 2, 
      plot_designs = F,
      method = "finite", n_cox_points = 1000,
      verbose = 0,
      max_it = 2,
      seed = 100)
  },
  
  e_R = {
    res_alg_order_R = coord_ex_mixt_R(
      20, q = 15, 
      order = 2, 
      plot_designs = F,
      method = "finite", n_cox_points = 1000,
      verbose = 0,
      max_it = 2,
      seed = 100)
  },
  times = 1
)



# Takes like a minute
# 8x speedup
time_res_2 = microbenchmark(
  e_C = {
    res_alg_order_C = coord_ex_mixt(
      10, q = 5, 
      order = 2, 
      plot_designs = F,
      method = "finite", n_cox_points = 1000,
      verbose = 0,
      max_it = 2,
      seed = 100)
  },
  
  e_R = {
    res_alg_order_R = coord_ex_mixt_R(
      10, q = 5, 
      order = 2, 
      plot_designs = F,
      method = "finite", n_cox_points = 1000,
      verbose = 0,
      max_it = 2,
      seed = 100)
  },
  times = 2
)



# takes like 4 minutes
# 3.7x speedup
time_res_3 = microbenchmark(
  e_C = {
    res_alg_order_C = coord_ex_mixt(
      100, q = 5, 
      order = 2, 
      plot_designs = F,
      method = "finite", n_cox_points = 1000,
      verbose = 0,
      max_it = 2,
      seed = 100)
  },
  
  e_R = {
    res_alg_order_R = coord_ex_mixt_R(
      100, q = 5, 
      order = 2, 
      plot_designs = F,
      method = "finite", n_cox_points = 1000,
      verbose = 0,
      max_it = 2,
      seed = 100)
  },
  times = 1
)








init_time = Sys.time()
# results_coord_ex_all = lapply(c(10, 100, 1000, 10000), function(n_runs){
results_coord_ex_small = lapply(c(10, 100), function(n_runs){
  # lapply(seq(5, 20, by = 5), function(q){
  lapply(c(3, 5), function(q){
    lapply(c(100, 1000), function(n_cox_points){
      lapply(1:3, function(order){
        seed = as.numeric(paste0(n_cox_points, order))
        
        cat(
          "\n\n\n\n",
          "n_runs = ", n_runs,
          ", q = ", q,
          ", n_cox_points = ", n_cox_points,
          ", order = ", order,
          ", seed = ", seed, "\n\t",
          as.character(Sys.time()),
          "\n",
          sep = ""
        )
        
        time_res = microbenchmark(
          e_C = {
            res_alg_order_C = coord_ex_mixt(
              n_runs, q = q, order = order, 
              plot_designs = F,
              method = "finite", n_cox_points = n_cox_points,
              verbose = 0,
              max_it = 2,
              seed = seed)
          },
          
          e_R = {
            res_alg_order_R = coord_ex_mixt_R(
              n_runs, q = q, 
              order = order, 
              plot_designs = F,
              method = "finite", n_cox_points = n_cox_points,
              verbose = 0,
              max_it = 2,
              seed = seed)
          },
          times = 2
        )
        
        out_df = time_res %>% 
          as.data.frame() %>% 
          group_by(expr) %>% 
          summarize(
            mean = mean(time),
            median = median(time),
            p25 = as.numeric(quantile(time, 0.25)),
            p75 = as.numeric(quantile(time, 0.75))) %>% 
          mutate(n_runs = n_runs,
                 q = q,
                 n_cox_points = n_cox_points,
                 order = order)
        
        return(out_df)
        
      }) %>% 
        bind_rows()
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

end_time = Sys.time()

end_time - init_time


saveRDS(results_coord_ex_small, "results_coord_ex_small.rds")





df_plot_small = results_coord_ex_small %>% 
  select(median, expr, n_runs, q, n_cox_points, order) %>% 
  pivot_wider(values_from = c(median), 
              names_from = expr) %>% 
  mutate(median_speedup = e_R/e_C,
         q = as.factor(q),
         n_runs = as.factor(n_runs),
         n_cox_points = as.factor(n_cox_points),
         order = as.factor(order)
         ) %>% 
  select(n_runs, 
         q, 
         median_speedup,
         n_cox_points,
         order) %>% 
  inner_join(
    results_coord_ex_small %>% 
      select(mean, expr, n_runs, q, n_cox_points, order) %>% 
      pivot_wider(values_from = c(mean), 
                  names_from = expr) %>% 
      mutate(mean_speedup = e_R/e_C,
             q = as.factor(q),
             n_runs = as.factor(n_runs),
             n_cox_points = as.factor(n_cox_points),
             order = as.factor(order)
      ) %>% 
      select(n_runs, 
             q, 
             mean_speedup,
             n_cox_points,
             order)
  )




df_plot_small %>%
  ggplot(aes(n_cox_points, median_speedup, color = order, group = order)) +
  geom_point() +
  geom_line() +
  facet_grid(n_runs~q, scales = "free")







init_time = Sys.time()
# results_coord_ex_all = lapply(c(10, 100, 1000, 10000), function(n_runs){
results_coord_ex_all = lapply(c(10, 100), function(n_runs){
  # lapply(seq(5, 20, by = 5), function(q){
  lapply(seq(5, 15, by = 5), function(q){
    lapply(c(10, 100, 1000), function(n_cox_points){
      lapply(1:3, function(order){
        seed = as.numeric(paste0(n_cox_points, order))
        
        cat(
          "\n\n\n\n",
          "n_runs = ", n_runs,
          ", q = ", q,
          ", n_cox_points = ", n_cox_points,
          ", order = ", order,
          ", seed = ", seed, "\n\t",
          as.character(Sys.time()),
          "\n",
          sep = ""
        )
        
        time_res = microbenchmark(
          e_C = {
            res_alg_order_C = coord_ex_mixt(
              n_runs, q = q, order = order, 
              plot_designs = F,
              method = "finite", n_cox_points = n_cox_points,
              verbose = 0,
              max_it = 2,
              seed = seed)
          },
          
          e_R = {
            res_alg_order_R = coord_ex_mixt_R(
              n_runs, q = q, 
              order = order, 
              plot_designs = F,
              method = "finite", n_cox_points = n_cox_points,
              verbose = 0,
              max_it = 2,
              seed = seed)
          },
          times = 2
        )
        
        out_df = time_res %>% 
          as.data.frame() %>% 
          group_by(expr) %>% 
          summarize(
            mean = mean(time),
            median = median(time),
            p25 = as.numeric(quantile(time, 0.25)),
            p75 = as.numeric(quantile(time, 0.75))) %>% 
          mutate(n_runs = n_runs,
                 q = q,
                 n_cox_points = n_cox_points,
                 order = order)
        
        return(out_df)
        
      }) %>% 
        bind_rows()
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()


end_time = Sys.time()

end_time - init_time


saveRDS(results_coord_ex_all, "results_coord_ex_all.rds")
saveRDS(end_time - init_time, "time.rds")

# It started at "2019-12-18 23:19:22 CET"
# n_runs = 100, q = 5, n_cox_points = 100, order = 1, seed = 1001
# 2019-12-19 21:39:37
# 
# n_runs = 100, q = 5, n_cox_points = 100, order = 2, seed = 1002
# 2019-12-19 21:39:58
#  Error in compute_cox_direction(X[k, ], i, n_cox_points) : 
#   Error while computing Cox direction. Value out of bounds.
# Cox direction computed:
# 	c(0, 0, 0, -0.127551020408163, 0) 
#
# Ended at "2019-12-19 21:40:03 CET" because of the previous error.

