# function that loops through a moving window of X days of tracking data
# computes estimated HR size and CIs for each window
# plots ML point estimate of HR size and 95% CIs
# written by Stefano Mezzini based on code by Mark Bidwell and Chris Fleming

library('ctmm')  # for continuous movement modeling
library('dplyr') # for data wrangling
library('purrr') # for functional programming

window_hr <- function(tel, window, dt, projection, full_ud = NULL,
                      fig_path = NULL, rds_path = NULL, cores = 1,
                      weights = FALSE) {
  
  if(! is.null(full_ud) & class(full_ud) == 'UD') {
    HR_0 <- summary(full_ud, units = FALSE)$CI['area (square meters)',
                                                 'est'] / 1e6
  } else {
    HR_0 <- NA_real_
  }
  
  # moving window created at beginning of t; slides forward as far as possible:
  # |---|.... --> .|---|... --> ..|---|.. --> ...|---|. --> ....|---|
  
  times <- seq(min(tel$t), max(tel$t) - window, by = dt)
  N <- length(times)
  
  # to extract home ranges later
  extract_hr <- function(a, par, l.ud) {
    # convert HR to km^2
    summary(a, units = FALSE, level.UD = l.ud)$CI['area (square meters)', ][par] / 1e6
  }
  
  out <-
    tibble(
      # add start and end times
      t_start = times, # left bound
      t_end = t_start + window, # right bound
      # subset times within window (t to t + win; can't use filter() on telem.)
      dataset = map2(t_start, t_end,
                     function(t_1, t_2) tel[tel$t >= t_1 & tel$t <= t_2, ]),
      models =
        imap(dataset,
             \(d, i) {
               
               if(nrow(d) > 1) {
                 cat('Analyzing dataset ', i, ' of ', N, '.\n', sep = '')
                 tibble(
                   # find initial guesses for models
                   guess = ctmm.guess(data = d, interactive = FALSE) %>% list(),
                   # select best model based on subset of tel
                   model = ctmm.select(data = d, CTMM = guess[[1]],
                                       cores = cores) %>%
                     list(),
                   # estimate autocorrelated kernel density estimate
                   akde = akde(data = d, CTMM = model[[1]], weights = weights) %>% list(),
                   # find home range estimate
                   hr_est_50 = extract_hr(a = akde[[1]], par='est', l.ud=0.50),
                   hr_lwr_50 = extract_hr(a = akde[[1]], par='low', l.ud=0.50),
                   hr_upr_50 = extract_hr(a = akde[[1]], par='high', l.ud=0.50),
                   hr_est_95 = extract_hr(a = akde[[1]], par='est', l.ud=0.95),
                   hr_lwr_95 = extract_hr(a = akde[[1]], par='low', l.ud=0.95),
                   hr_upr_95 = extract_hr(a = akde[[1]], par='high', l.ud=0.95))
               } else {
                 tibble(
                   # find initial guesses for models
                   guess = list('Insufficient data.'),
                   # select best model based on subset of tel
                   model = list('Insufficient data.'),
                   # estimate autocorrelated kernel density estimate
                   akde = list('Insufficient data.'),
                   # find home range estimate
                   hr_est_50 = NA_real_,
                   hr_lwr_50 = NA_real_,
                   hr_upr_50 = NA_real_,
                   hr_est_95 = NA_real_,
                   hr_lwr_95 = NA_real_,
                   hr_upr_95 = NA_real_)
               } # close else
             })) %>% # close function for imap()
    tidyr::unnest(models) %>%
    # mutate(ess = map_dbl(model, \(.m) summary(.m)$DOF['area'])) %>%
    mutate(t_center = (t_start + t_end) / 2,
           posixct = as.POSIXct(t_center, origin = '1970-01-01',
                                tz = tel@info$timezone),
           date = as.Date(posixct))
  
  if(! is.null(rds_path)) {
    saveRDS(out, file.path(rds_path,
                           paste0(tel@info['identity'],
                                  '-window-', window / (1 %#% 'day'), '-days',
                                  '-dt-', dt / (1 %#% 'day'), '-days.rds')))
  }
  
  # plot results
  library('ggplot2') # for fancy plots
  theme_set(theme_bw() + theme(legend.position = 'none'))
  
  # tracking data
  plt_a <-
    ggplot(tel) +
    coord_equal() +
    geom_point(aes(longitude, latitude, color = timestamp)) +
    geom_path(aes(longitude, latitude), alpha = 0.1) +
    scale_color_viridis_c('Time') +
    labs(x = '', y = NULL)
  
  plt_b <-
    ggplot(out) +
    
    # # 95% CIs for home ranges
    geom_ribbon(aes(date, ymin = hr_lwr_50, ymax = hr_upr_50), alpha = 0.3) +
    geom_ribbon(aes(date, ymin = hr_lwr_95, ymax = hr_upr_95), alpha = 0.3) +
    
    # # core home range
    # geom_line(aes(date, hr_est_50), linewidth = 1.25) +
    # geom_line(aes(date, hr_est_50, color = posixct)) +
    
    # 95% home range
    geom_line(aes(date, hr_est_95), linewidth = 1.25) +
    geom_line(aes(date, hr_est_95, color = posixct)) +
    
    # home range from model fit to full dataset
    geom_hline(yintercept = HR_0, color = 'darkorange') +
    
    scale_x_date(NULL, date_labels = '%b %Y') +
    scale_color_viridis_c() +
    labs(y = expression(Home~range~(km^2)))
  
  plt <- cowplot::plot_grid(plt_a, plt_b, labels = c('a.', 'b.'), nrow = 1,
                            align = 'hv')
  
  if (! is.null(fig_path)) {
    # Save figure as a png using the animal's name
    file.path(fig_path, paste0(tel@info['identity'],
                               '-window-', window / (1 %#% 'day'),
                               '-days-dt-', dt / (1 %#% 'day'),
                               '-days.png')) %>%
      ggsave(plot = plt, units = 'in', width = 7, height = 3,
             dpi = 600, bg = 'white')
  } else {
    print(plt)
  }
}

if(FALSE) {
  # Test the function on the buffalo data
  data('buffalo')
  
  # Apply to a single individual
  test <- window_hr(tel = buffalo$Queen, # buffalo with smallest dataset
                    window = 7 %#% 'day', # size of the moving window
                    dt = 1 %#% 'day', # step size for the moving window
                    fig_path = NULL) # can specify where to save the figure
  
  # apply to all buffaloes in the dataset
  results <- lapply(buffalo,
                    window_hr,
                    # arguments to pass to window_hr
                    window = 30 %#% 'day',
                    dt = 2 %#% 'day')
}
