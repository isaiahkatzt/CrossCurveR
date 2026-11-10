test_that("withdrawal sampling selects the sole eligible survivor", {
  state <- tontine_initialize_state(tontine_create_member_pool(7, initial_position = 100))
  for (survivor in c(1L, 7L)) {
    alive <- state
    alive$members$alive <- alive$members$member_id == survivor
    for (eligible in list(NULL, survivor)) {
      set.seed(123)
      ids <- tontine_sample_withdrawal_ids(alive, 100, eligible)
      expect_identical(ids, survivor)
      set.seed(123)
      result <- tontine_sample_withdrawal_ids_from_integrated_intensity(
        alive, 100, eligible, warn_on_cap = FALSE)
      expect_identical(result$withdrawal_member_ids, survivor)
      expect_equal(result$withdrawal_count, 1)
      expect_true(result$withdrawal_count_capped)
      applied <- tontine_apply_withdrawals(alive, 1, result$withdrawal_member_ids)
      expect_identical(applied$event$member_id, survivor)
    }
  }
  # Explicit post-death eligibility may be narrower than the current state.
  set.seed(123)
  expect_identical(tontine_sample_withdrawal_ids(state, 100, 7L), 7L)
  set.seed(123)
  expect_identical(tontine_sample_withdrawal_ids_from_integrated_intensity(
    state, 100, 7L, warn_on_cap = FALSE)$withdrawal_member_ids, 7L)
})

test_that("withdrawal sampling handles empty eligibility and zero intensity", {
  state <- tontine_initialize_state(tontine_create_member_pool(7))
  for (eligible in list(integer(), 7L, c(2L, 7L))) {
    intensity <- if (length(eligible)) 0 else 100
    expect_identical(tontine_sample_withdrawal_ids(state, intensity, eligible), integer())
    result <- tontine_sample_withdrawal_ids_from_integrated_intensity(state, intensity, eligible)
    expect_identical(result$withdrawal_member_ids, integer())
    expect_equal(result$withdrawal_count, 0)
    expect_false(result$withdrawal_count_capped)
  }
})
